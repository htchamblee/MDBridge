// MP4 encoding with the FFmpeg libraries; see encoder.hpp. Port of MDATK's FFmpegWriter
// (render/video.py, owner's code) from a piped `ffmpeg` process to libavcodec / libavformat.
#include "video/encoder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <mutex>
#include <string>
#include <system_error>
#include <utility>

#include "mdbridge/core/paths.hpp"
#include "video/ffmpeg_support.hpp"

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

#if MDBRIDGE_HAVE_FFMPEG
extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/avutil.h>
#include <libavutil/dict.h>
#include <libavutil/frame.h>
#include <libavutil/log.h>
#include <libavutil/opt.h>
#include <libavutil/pixdesc.h>
#include <libswscale/swscale.h>
}
#endif

namespace fs = std::filesystem;

namespace mdbridge::render::video {

namespace {

constexpr const char* kNotAvailable = "video encoding is not available in this build";

// SVT-AV1 speed / quality trade-off (0 slowest .. 13 fastest); 8 encodes 1080p animations at
// interactive speed with good quality.
constexpr int kSvtAv1Preset = 8;

int clampCrf(int crf) { return std::clamp(crf, 0, 51); }

long processId() {
#ifdef _WIN32
    return static_cast<long>(_getpid());
#else
    return static_cast<long>(getpid());
#endif
}

}  // namespace

// ---- codec choice ----

int svtAv1Crf(int crf) {
    const int c = clampCrf(crf);
    return (c * 63 + 25) / 51;  // round(c * 63 / 51), exact in integers
}

int mpeg4Qscale(int crf) {
    const int c = clampCrf(crf);
    if (c <= 20) {
        return 1 + (c + 5) / 10;  // round(1 + c / 10): 0 -> 1, 10 -> 2, 20 -> 3
    }
    return 3 + ((c - 20) * 28 + 15) / 31;  // 21 -> 4, 51 -> 31
}

CodecChoice chooseCodec(std::string_view requested, int crf,
                        const std::function<bool(std::string_view)>& available) {
    CodecChoice choice;
    choice.requested = requested.empty() ? std::string("libx264") : std::string(requested);
    const std::string& name = choice.requested;

    // the H.264 / H.265 names (GPL encoders or patent-encumbered formats) and the AV1 family name
    // mean "the default encoder": AV1 (owner decision D1)
    const bool gplEncoder = name == "libx264" || name == "libx265";
    const bool aliasOfAv1 = gplEncoder || name == "h264" || name == "hevc" || name == "av1";

    auto use = [&](std::string_view encoder) {
        choice.encoder = std::string(encoder);
        if (encoder == kAv1Encoder) {
            choice.qualityOption = "crf";
            choice.quality = svtAv1Crf(crf);
        } else if (encoder == kFallbackEncoder) {
            choice.qualityOption = "qscale";
            choice.quality = mpeg4Qscale(crf);
        }
    };

    if (aliasOfAv1 && available(kAv1Encoder)) {
        use(kAv1Encoder);
        if (name != "av1") {
            choice.message = "video: codec " + name + " is not bundled" + (gplEncoder ? " (GPL)" : "") +
                             "; using " + std::string(kAv1Encoder) + " (AV1)";
        }
        return choice;
    }
    if (!aliasOfAv1 && available(name)) {
        use(name);
        return choice;
    }
    for (std::string_view fallback : {kAv1Encoder, kFallbackEncoder}) {
        if (available(fallback)) {
            use(fallback);
            choice.message = "video: encoder " + name + " is not available; using " + std::string(fallback);
            return choice;
        }
    }
    throw VideoError("no video encoder is available (tried " + name + ", " + std::string(kAv1Encoder) +
                     " and " + std::string(kFallbackEncoder) + ")");
}

fs::path temporaryVideoPath(const fs::path& path) {
    const fs::path absolute = fs::absolute(path);
    return absolute.parent_path() /
           core::pathFromUtf8(".MDATK_tmp_" + std::to_string(processId()) + "_" +
                              core::pathToUtf8(absolute.filename()));
}

#if MDBRIDGE_HAVE_FFMPEG

// ---- FFmpeg support ----

namespace detail {

namespace {

std::mutex g_logMutex;
std::deque<std::string> g_log;  // the last few error messages

void logCallback(void* avcl, int level, const char* fmt, va_list args) {
    if (level > AV_LOG_ERROR) {
        return;
    }
    std::array<char, 1024> line{};
    int printPrefix = 1;
    av_log_format_line2(avcl, level, fmt, args, line.data(), static_cast<int>(line.size()), &printPrefix);
    std::string text(line.data());
    // "[libsvtav1 @ 0x55d0c1a2b3c0] ..." -> "[libsvtav1] ...": no addresses in messages
    const std::size_t at = text.find(" @ 0x");
    const std::size_t close = at == std::string::npos ? at : text.find(']', at);
    if (close != std::string::npos) {
        text.erase(at, close - at);
    }
    while (!text.empty() && (text.back() == '\n' || text.back() == '\r' || text.back() == ' ')) {
        text.pop_back();
    }
    if (text.empty()) {
        return;
    }
    const std::lock_guard<std::mutex> lock(g_logMutex);
    g_log.push_back(std::move(text));
    while (g_log.size() > 8) {
        g_log.pop_front();
    }
}

}  // namespace

void initFfmpeg() {
    static std::once_flag once;
    std::call_once(once, [] {
        av_log_set_callback(logCallback);
        // SVT-AV1 prints its settings, and on abort() "deinit called without sending EOS", to
        // stderr; its real errors reach FFmpeg's log (and so the VideoError messages) anyway
        if (std::getenv("SVT_LOG") == nullptr) {
#ifdef _WIN32
            _putenv_s("SVT_LOG", "0");
#else
            setenv("SVT_LOG", "0", 0);
#endif
        }
    });
}

std::string takeFfmpegLog() {
    const std::lock_guard<std::mutex> lock(g_logMutex);
    std::string text;
    for (const std::string& line : g_log) {
        text += (text.empty() ? "" : "; ") + line;
    }
    g_log.clear();
    return text;
}

std::string ffmpegError(int code) {
    std::array<char, 256> buffer{};
    if (av_strerror(code, buffer.data(), buffer.size()) < 0) {
        std::snprintf(buffer.data(), buffer.size(), "error %d", code);
    }
    std::string text(buffer.data());
    const std::string log = takeFfmpegLog();
    if (!log.empty()) {
        text += " (" + log + ")";
    }
    return text;
}

}  // namespace detail

namespace {

std::string versionText(unsigned version) {
    return std::to_string(version >> 16) + "." + std::to_string((version >> 8) & 0xFFu) + "." +
           std::to_string(version & 0xFFu);
}

bool hasToken(const std::string& text, std::string_view token) {
    std::size_t pos = 0;
    while (pos < text.size()) {
        const std::size_t end = std::min(text.find(' ', pos), text.size());
        if (std::string_view(text).substr(pos, end - pos) == token) {
            return true;
        }
        pos = end + 1;
    }
    return false;
}

bool encoderExists(std::string_view name) {
    return avcodec_find_encoder_by_name(std::string(name).c_str()) != nullptr;
}

}  // namespace

bool videoEncodingAvailable() { return true; }

FfmpegBuildInfo ffmpegBuildInfo() {
    detail::initFfmpeg();
    FfmpegBuildInfo info;
    info.version = av_version_info();  // "n9.0.2" from a git checkout of the tag, "9.0.2" from a tarball
    if (info.version.size() > 1 && info.version[0] == 'n' && info.version[1] >= '0' && info.version[1] <= '9') {
        info.version.erase(0, 1);
    }
    info.configuration = avcodec_configuration();
    info.license = avcodec_license();
    info.libraries = "libavcodec " + versionText(avcodec_version()) + ", libavformat " +
                     versionText(avformat_version()) + ", libavutil " + versionText(avutil_version()) +
                     ", libswscale " + versionText(swscale_version());
    const bool lgplLicences = std::string_view(avcodec_license()).rfind("LGPL", 0) == 0 &&
                              std::string_view(avformat_license()).rfind("LGPL", 0) == 0 &&
                              std::string_view(avutil_license()).rfind("LGPL", 0) == 0 &&
                              std::string_view(swscale_license()).rfind("LGPL", 0) == 0;
    const std::string formatConfiguration = avformat_configuration();
    info.lgplOnly = lgplLicences && !hasToken(info.configuration, "--enable-gpl") &&
                    !hasToken(info.configuration, "--enable-nonfree") &&
                    !hasToken(formatConfiguration, "--enable-gpl") &&
                    !hasToken(formatConfiguration, "--enable-nonfree");
    info.hasSvtAv1 = encoderExists(kAv1Encoder);
    info.hasMpeg4 = encoderExists(kFallbackEncoder);
    return info;
}

std::string ffmpegVersionLine() {
    const FfmpegBuildInfo info = ffmpegBuildInfo();
    return "FFmpeg " + info.version + " libraries (" + info.license + "): " + info.libraries;
}

// ---- the encoder ----

struct VideoEncoder::Impl {
    fs::path path;
    fs::path temporary;
    EncoderOptions options;
    std::function<void(const std::string&)> warn;
    CodecChoice codec;
    AVRational rate{0, 1};
    AVFormatContext* format = nullptr;
    AVCodecContext* encoder = nullptr;
    AVStream* stream = nullptr;
    SwsContext* scaler = nullptr;
    AVFrame* frame = nullptr;
    AVPacket* packet = nullptr;
    int frames = 0;
    bool fileOpen = false;
    bool done = false;    // finish() succeeded
    bool failed = false;  // an error or abort(): nothing more can be written

    ~Impl() { release(); }

    void release() noexcept {
        sws_freeContext(scaler);
        scaler = nullptr;
        av_frame_free(&frame);
        av_packet_free(&packet);
        avcodec_free_context(&encoder);
        if (format != nullptr) {
            if (fileOpen) {
                avio_closep(&format->pb);
                fileOpen = false;
            }
            avformat_free_context(format);
            format = nullptr;
        }
        stream = nullptr;
    }

    [[noreturn]] void fail(const std::string& message) {
        failed = true;
        removeTemporary();
        throw VideoError(message);
    }

    void removeTemporary() noexcept {
        release();
        std::error_code ec;
        fs::remove(temporary, ec);
    }

    std::string sizeText() const {
        return std::to_string(options.width) + "x" + std::to_string(options.height);
    }

    // Opens the encoder `name` with the chosen quality; returns 0 or an FFmpeg error code.
    int openEncoder(const std::string& name) {
        avcodec_free_context(&encoder);
        const AVCodec* codecImpl = avcodec_find_encoder_by_name(name.c_str());
        if (codecImpl == nullptr) {
            return AVERROR_ENCODER_NOT_FOUND;
        }
        encoder = avcodec_alloc_context3(codecImpl);
        if (encoder == nullptr) {
            return AVERROR(ENOMEM);
        }
        encoder->width = options.width;
        encoder->height = options.height;
        encoder->pix_fmt = AV_PIX_FMT_YUV420P;
        encoder->framerate = rate;
        encoder->time_base = av_inv_q(rate);
        encoder->thread_count = std::max(0, options.threads);
        // sRGB frames converted with the BT.709 matrix to limited range, and tagged so: players
        // then show the colours as drawn (an untagged stream is read as BT.709 at HD sizes but
        // FFmpeg's default conversion is BT.601, which shifts hues)
        encoder->colorspace = AVCOL_SPC_BT709;
        encoder->color_primaries = AVCOL_PRI_BT709;
        encoder->color_trc = AVCOL_TRC_BT709;
        encoder->color_range = AVCOL_RANGE_MPEG;
        if ((format->oformat->flags & AVFMT_GLOBALHEADER) != 0) {
            encoder->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
        }
        if (codec.qualityOption == "crf") {
            av_opt_set_int(encoder->priv_data, "crf", codec.quality, 0);
            av_opt_set_int(encoder->priv_data, "preset", kSvtAv1Preset, 0);
        } else if (codec.qualityOption == "qscale") {  // `-q:v N`
            encoder->flags |= AV_CODEC_FLAG_QSCALE;
            encoder->global_quality = codec.quality * FF_QP2LAMBDA;
        }
        return avcodec_open2(encoder, codecImpl, nullptr);
    }

    void open() {
        if (options.width < 2 || options.height < 2 || options.width % 2 != 0 || options.height % 2 != 0) {
            throw VideoError("video size " + sizeText() + " is not even (yuv420p needs even sizes of at least 2x2)");
        }
        if (!(options.fps > 0.0) || !std::isfinite(options.fps)) {
            throw VideoError("video frame rate " + std::to_string(options.fps) + " is not positive");
        }
        detail::initFfmpeg();
        const FfmpegBuildInfo info = ffmpegBuildInfo();
        if (!info.lgplOnly) {
            throw VideoError("the FFmpeg libraries in use are not an LGPL-only build (" + info.license +
                             "; configuration: " + info.configuration +
                             "); MDBridge encodes video only with LGPL builds of FFmpeg");
        }
        codec = chooseCodec(options.codec, options.crf, encoderExists);
        if (!codec.message.empty() && warn) {
            warn(codec.message);
        }
        // the frame rate as `ffmpeg -r <fps>` reads it (av_parse_video_rate)
        rate = av_d2q(options.fps, 1001000);
        if (rate.num <= 0 || rate.den <= 0) {
            throw VideoError("video frame rate " + std::to_string(options.fps) + " is not usable");
        }

        temporary = temporaryVideoPath(path);
        const std::string temporaryUtf8 = core::pathToUtf8(temporary);
        int err = avformat_alloc_output_context2(&format, nullptr, nullptr, temporaryUtf8.c_str());
        if (err < 0 || format == nullptr) {
            detail::takeFfmpegLog();  // "unable to find a suitable output format": use mp4
            err = avformat_alloc_output_context2(&format, nullptr, "mp4", temporaryUtf8.c_str());
        }
        if (err < 0 || format == nullptr) {
            fail("video: no MP4 writer in the FFmpeg libraries: " + detail::ffmpegError(err));
        }

        err = openEncoder(codec.encoder);
        if (err < 0 && codec.encoder != kFallbackEncoder && encoderExists(kFallbackEncoder)) {
            const std::string reason = detail::ffmpegError(err);
            const std::string tried = codec.encoder;
            codec.encoder = std::string(kFallbackEncoder);
            codec.qualityOption = "qscale";
            codec.quality = mpeg4Qscale(options.crf);
            const std::string message = "video: " + tried + " cannot encode " + sizeText() + " (" + reason +
                                        "); using " + codec.encoder;
            codec.message = codec.message.empty() ? message : codec.message + "; " + message;
            if (warn) {
                warn(message);
            }
            err = openEncoder(codec.encoder);
        }
        if (err < 0) {
            fail("video: cannot open the encoder " + codec.encoder + " for " + sizeText() + " at " +
                 std::to_string(rate.num) + "/" + std::to_string(rate.den) + " fps: " + detail::ffmpegError(err));
        }

        stream = avformat_new_stream(format, nullptr);
        if (stream == nullptr) {
            fail("video: cannot add a video stream: out of memory");
        }
        stream->time_base = encoder->time_base;  // the muxer may refine it (mp4: x 2^k >= 10000)
        stream->avg_frame_rate = rate;
        err = avcodec_parameters_from_context(stream->codecpar, encoder);
        if (err < 0) {
            fail("video: " + detail::ffmpegError(err));
        }

        if ((format->oformat->flags & AVFMT_NOFILE) == 0) {
            err = avio_open(&format->pb, temporaryUtf8.c_str(), AVIO_FLAG_WRITE);
            if (err < 0) {
                fail("video: cannot write " + temporaryUtf8 + ": " + detail::ffmpegError(err));
            }
            fileOpen = true;
        }
        AVDictionary* muxerOptions = nullptr;
        av_dict_set(&muxerOptions, "movflags", "+faststart", 0);
        err = avformat_write_header(format, &muxerOptions);
        av_dict_free(&muxerOptions);
        if (err < 0) {
            fail("video: cannot start " + temporaryUtf8 + ": " + detail::ffmpegError(err));
        }

        scaler = sws_getContext(options.width, options.height, AV_PIX_FMT_RGBA, options.width, options.height,
                                AV_PIX_FMT_YUV420P, SWS_BICUBIC | SWS_ACCURATE_RND, nullptr, nullptr, nullptr);
        if (scaler == nullptr) {
            fail("video: cannot convert RGBA to yuv420p at " + sizeText() + ": " + detail::takeFfmpegLog());
        }
        // source: full-range RGB; destination: BT.709, limited range
        sws_setColorspaceDetails(scaler, sws_getCoefficients(SWS_CS_ITU709), 1,
                                 sws_getCoefficients(SWS_CS_ITU709), 0, 0, 1 << 16, 1 << 16);

        frame = av_frame_alloc();
        packet = av_packet_alloc();
        if (frame == nullptr || packet == nullptr) {
            fail("video: out of memory");
        }
        frame->format = AV_PIX_FMT_YUV420P;
        frame->width = options.width;
        frame->height = options.height;
        frame->colorspace = encoder->colorspace;
        frame->color_primaries = encoder->color_primaries;
        frame->color_trc = encoder->color_trc;
        frame->color_range = encoder->color_range;
        err = av_frame_get_buffer(frame, 0);
        if (err < 0) {
            fail("video: " + detail::ffmpegError(err));
        }
    }

    // Sends `input` (nullptr: flush) and writes every packet the encoder returns.
    void encode(const AVFrame* input) {
        int err = avcodec_send_frame(encoder, input);
        if (err < 0) {
            fail("video: the encoder " + codec.encoder + " failed at frame " + std::to_string(frames + 1) + ": " +
                 detail::ffmpegError(err));
        }
        for (;;) {
            err = avcodec_receive_packet(encoder, packet);
            if (err == AVERROR(EAGAIN) || err == AVERROR_EOF) {
                return;
            }
            if (err < 0) {
                fail("video: the encoder " + codec.encoder + " failed: " + detail::ffmpegError(err));
            }
            packet->stream_index = stream->index;
            if (packet->duration <= 0) {
                packet->duration = 1;  // one frame in the encoder's time base 1 / fps
            }
            av_packet_rescale_ts(packet, encoder->time_base, stream->time_base);
            err = av_interleaved_write_frame(format, packet);  // takes the packet's data
            if (err < 0) {
                fail("video: cannot write " + core::pathToUtf8(temporary) + ": " + detail::ffmpegError(err));
            }
        }
    }
};

VideoEncoder::VideoEncoder(const fs::path& path, const EncoderOptions& options,
                           std::function<void(const std::string&)> warn)
    : m(std::make_unique<Impl>()) {
    m->path = path;
    m->options = options;
    m->warn = std::move(warn);
    m->open();
}

VideoEncoder::~VideoEncoder() {
    if (m && !m->done) {
        abort();
    }
}

void VideoEncoder::writeFrame(const std::uint8_t* rgba, std::size_t strideBytes) {
    if (m->done || m->failed) {
        throw VideoError("video: the encoder is closed");
    }
    if (rgba == nullptr) {
        throw VideoError("video: no frame data");
    }
    const std::size_t stride = strideBytes == 0 ? static_cast<std::size_t>(m->options.width) * 4 : strideBytes;
    if (stride < static_cast<std::size_t>(m->options.width) * 4 || stride > 0x7FFFFFFF) {
        throw VideoError("video: frame row stride " + std::to_string(stride) + " is too small for width " +
                         std::to_string(m->options.width));
    }
    const int err = av_frame_make_writable(m->frame);
    if (err < 0) {
        m->fail("video: " + detail::ffmpegError(err));
    }
    const std::uint8_t* const sourceData[4] = {rgba, nullptr, nullptr, nullptr};
    const int sourceStride[4] = {static_cast<int>(stride), 0, 0, 0};
    sws_scale(m->scaler, sourceData, sourceStride, 0, m->options.height, m->frame->data, m->frame->linesize);
    m->frame->pts = m->frames;
    m->frame->duration = 1;
    if (m->codec.qualityOption == "qscale") {
        m->frame->quality = m->encoder->global_quality;  // as ffmpeg does for `-q:v`
    }
    m->encode(m->frame);
    ++m->frames;
}

fs::path VideoEncoder::finish() {
    if (m->done) {
        return m->path;
    }
    if (m->failed) {
        throw VideoError("video: the encoder is closed");
    }
    if (m->frames == 0) {
        m->fail("video: no frames to write");
    }
    m->encode(nullptr);
    int err = av_write_trailer(m->format);  // faststart: moves moov in front of mdat here
    if (err < 0) {
        m->fail("video: cannot finish " + core::pathToUtf8(m->temporary) + ": " + detail::ffmpegError(err));
    }
    if (m->fileOpen) {
        err = avio_closep(&m->format->pb);
        m->fileOpen = false;
        if (err < 0) {
            m->fail("video: cannot write " + core::pathToUtf8(m->temporary) + ": " + detail::ffmpegError(err));
        }
    }
    m->release();
    std::error_code ec;
    fs::rename(m->temporary, m->path, ec);  // replaces an existing file (os.replace)
    if (ec) {
        m->fail("video: cannot rename " + core::pathToUtf8(m->temporary) + " to " + core::pathToUtf8(m->path) +
                ": " + ec.message());
    }
    m->done = true;
    return m->path;
}

void VideoEncoder::abort() noexcept {
    if (m->done) {
        return;
    }
    m->failed = true;
    m->removeTemporary();
}

int VideoEncoder::framesWritten() const { return m->frames; }
const CodecChoice& VideoEncoder::codec() const { return m->codec; }
const fs::path& VideoEncoder::temporaryPath() const { return m->temporary; }
const fs::path& VideoEncoder::path() const { return m->path; }

#else  // no FFmpeg libraries in this build

bool videoEncodingAvailable() { return false; }
FfmpegBuildInfo ffmpegBuildInfo() { throw VideoError(kNotAvailable); }
std::string ffmpegVersionLine() { return "not available"; }

struct VideoEncoder::Impl {
    fs::path path;
    fs::path temporary;
    CodecChoice codec;
};

VideoEncoder::VideoEncoder(const fs::path&, const EncoderOptions&, std::function<void(const std::string&)>) {
    throw VideoError(kNotAvailable);
}
VideoEncoder::~VideoEncoder() = default;
void VideoEncoder::writeFrame(const std::uint8_t*, std::size_t) { throw VideoError(kNotAvailable); }
fs::path VideoEncoder::finish() { throw VideoError(kNotAvailable); }
void VideoEncoder::abort() noexcept {}
int VideoEncoder::framesWritten() const { return 0; }
const CodecChoice& VideoEncoder::codec() const { return m->codec; }
const fs::path& VideoEncoder::temporaryPath() const { return m->temporary; }
const fs::path& VideoEncoder::path() const { return m->path; }

#endif

}  // namespace mdbridge::render::video
