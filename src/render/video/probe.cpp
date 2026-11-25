// Demux-based self-check of a written video; see probe.hpp.
#include "video/probe.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "mdbridge/core/paths.hpp"
#include "video/encoder.hpp"
#include "video/ffmpeg_support.hpp"

#if MDBRIDGE_HAVE_FFMPEG
extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/pixdesc.h>
#include <libswscale/swscale.h>
}
#endif

namespace fs = std::filesystem;

namespace mdbridge::render::video {

namespace {

// Top-level ISO BMFF boxes in file order: is `moov` before the first `mdat`?
[[maybe_unused]] bool moovBeforeMdat(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }
    in.seekg(0, std::ios::end);
    const auto fileSize = static_cast<std::uint64_t>(in.tellg());
    std::uint64_t offset = 0;
    while (offset + 8 <= fileSize) {
        in.seekg(static_cast<std::streamoff>(offset));
        std::array<unsigned char, 16> header{};
        in.read(reinterpret_cast<char*>(header.data()), 8);
        if (!in) {
            return false;
        }
        std::uint64_t size = (std::uint64_t{header[0]} << 24) | (std::uint64_t{header[1]} << 16) |
                             (std::uint64_t{header[2]} << 8) | std::uint64_t{header[3]};
        const std::string type(reinterpret_cast<const char*>(header.data() + 4), 4);
        if (size == 1) {  // 64-bit size follows
            in.read(reinterpret_cast<char*>(header.data() + 8), 8);
            if (!in) {
                return false;
            }
            size = 0;
            for (int k = 8; k < 16; ++k) {
                size = (size << 8) | header[static_cast<std::size_t>(k)];
            }
        } else if (size == 0) {  // to the end of the file
            size = fileSize - offset;
        }
        if (type == "moov") {
            return true;
        }
        if (type == "mdat" || size < 8) {
            return false;
        }
        offset += size;
    }
    return false;
}

}  // namespace

#if MDBRIDGE_HAVE_FFMPEG

namespace {

constexpr std::int64_t kNoPts = std::numeric_limits<std::int64_t>::min();  // AV_NOPTS_VALUE

std::string nameOr(const char* name) { return name != nullptr ? std::string(name) : std::string("unknown"); }

// Pixel format from an AV1CodecConfigurationRecord (the MP4 `av1C` box, FFmpeg's extradata) when
// the libraries have no AV1 decoder to report it: bit depth, monochrome and chroma subsampling.
std::string av1PixelFormat(const AVCodecParameters* par) {
    if (par->extradata == nullptr || par->extradata_size < 4 || par->extradata[0] != 0x81) {
        return "unknown";
    }
    const unsigned flags = par->extradata[2];
    const bool high = (flags & 0x40u) != 0;
    const bool twelve = (flags & 0x20u) != 0;
    const bool mono = (flags & 0x10u) != 0;
    const bool subX = (flags & 0x08u) != 0;
    const bool subY = (flags & 0x04u) != 0;
    const int depth = !high ? 8 : (twelve ? 12 : 10);
    AVPixelFormat fmt = AV_PIX_FMT_NONE;
    if (mono) {
        fmt = depth == 8 ? AV_PIX_FMT_GRAY8 : (depth == 10 ? AV_PIX_FMT_GRAY10LE : AV_PIX_FMT_GRAY12LE);
    } else if (subX && subY) {
        fmt = depth == 8 ? AV_PIX_FMT_YUV420P : (depth == 10 ? AV_PIX_FMT_YUV420P10LE : AV_PIX_FMT_YUV420P12LE);
    } else if (subX) {
        fmt = depth == 8 ? AV_PIX_FMT_YUV422P : (depth == 10 ? AV_PIX_FMT_YUV422P10LE : AV_PIX_FMT_YUV422P12LE);
    } else if (!subY) {
        fmt = depth == 8 ? AV_PIX_FMT_YUV444P : (depth == 10 ? AV_PIX_FMT_YUV444P10LE : AV_PIX_FMT_YUV444P12LE);
    }
    return nameOr(av_get_pix_fmt_name(fmt));
}

struct InputCloser {
    AVFormatContext* format = nullptr;
    AVPacket* packet = nullptr;
    ~InputCloser() {
        av_packet_free(&packet);
        avformat_close_input(&format);
    }
};

}  // namespace

VideoProbe probeVideo(const fs::path& path) {
    detail::initFfmpeg();
    const std::string pathUtf8 = core::pathToUtf8(path);
    InputCloser input;
    int err = avformat_open_input(&input.format, pathUtf8.c_str(), nullptr, nullptr);
    if (err < 0) {
        throw VideoError("video: cannot open " + pathUtf8 + ": " + detail::ffmpegError(err));
    }
    err = avformat_find_stream_info(input.format, nullptr);
    if (err < 0) {
        throw VideoError("video: cannot read " + pathUtf8 + ": " + detail::ffmpegError(err));
    }
    VideoProbe probe;
    probe.format = nameOr(input.format->iformat->name);
    int videoIndex = -1;
    for (unsigned k = 0; k < input.format->nb_streams; ++k) {
        const AVCodecParameters* par = input.format->streams[k]->codecpar;
        if (par->codec_type == AVMEDIA_TYPE_VIDEO) {
            ++probe.videoStreams;
            if (videoIndex < 0) {
                videoIndex = static_cast<int>(k);
            }
        } else if (par->codec_type == AVMEDIA_TYPE_AUDIO) {
            ++probe.audioStreams;
        } else {
            ++probe.otherStreams;
        }
    }
    if (videoIndex >= 0) {
        const AVStream* stream = input.format->streams[videoIndex];
        const AVCodecParameters* par = stream->codecpar;
        probe.codec = nameOr(avcodec_get_name(par->codec_id));
        probe.pixelFormat = nameOr(av_get_pix_fmt_name(static_cast<AVPixelFormat>(par->format)));
        if (probe.pixelFormat == "unknown" && par->codec_id == AV_CODEC_ID_AV1) {
            probe.pixelFormat = av1PixelFormat(par);
        }
        probe.width = par->width;
        probe.height = par->height;
        probe.frameRateNum = stream->r_frame_rate.num;
        probe.frameRateDen = stream->r_frame_rate.den;
        probe.timeBaseNum = stream->time_base.num;
        probe.timeBaseDen = stream->time_base.den;
        probe.durationTs = stream->duration == kNoPts ? 0 : stream->duration;
        probe.durationSeconds = static_cast<double>(probe.durationTs) * av_q2d(stream->time_base);
        probe.nbFrames = stream->nb_frames;
        probe.colorSpace = nameOr(av_color_space_name(par->color_space));
        probe.colorRange = nameOr(av_color_range_name(par->color_range));
        probe.colorPrimaries = nameOr(av_color_primaries_name(par->color_primaries));
        probe.colorTransfer = nameOr(av_color_transfer_name(par->color_trc));

        input.packet = av_packet_alloc();
        if (input.packet == nullptr) {
            throw VideoError("video: out of memory");
        }
        bool first = true;
        while ((err = av_read_frame(input.format, input.packet)) >= 0) {
            if (input.packet->stream_index == videoIndex) {
                ++probe.packets;
                if ((input.packet->flags & AV_PKT_FLAG_KEY) != 0) {
                    ++probe.keyframes;
                }
                if (first) {
                    probe.firstPts = input.packet->pts == kNoPts ? 0 : input.packet->pts;
                    first = false;
                }
            }
            av_packet_unref(input.packet);
        }
        if (err != AVERROR_EOF) {
            throw VideoError("video: cannot read " + pathUtf8 + ": " + detail::ffmpegError(err));
        }
    }
    probe.moovBeforeMdat = moovBeforeMdat(path);
    return probe;
}

bool hasVideoDecoder(std::string_view codec) {
    const AVCodecDescriptor* desc = avcodec_descriptor_get_by_name(std::string(codec).c_str());
    return desc != nullptr && avcodec_find_decoder(desc->id) != nullptr;
}

namespace {

struct DecoderCloser {
    AVCodecContext* decoder = nullptr;
    AVFrame* frame = nullptr;
    SwsContext* scaler = nullptr;
    ~DecoderCloser() {
        sws_freeContext(scaler);
        av_frame_free(&frame);
        avcodec_free_context(&decoder);
    }
};

}  // namespace

std::vector<std::vector<std::uint8_t>> decodeVideoFrames(const fs::path& path) {
    detail::initFfmpeg();
    const std::string pathUtf8 = core::pathToUtf8(path);
    InputCloser input;
    int err = avformat_open_input(&input.format, pathUtf8.c_str(), nullptr, nullptr);
    if (err < 0) {
        throw VideoError("video: cannot open " + pathUtf8 + ": " + detail::ffmpegError(err));
    }
    err = avformat_find_stream_info(input.format, nullptr);
    if (err < 0) {
        throw VideoError("video: cannot read " + pathUtf8 + ": " + detail::ffmpegError(err));
    }
    const int index = av_find_best_stream(input.format, AVMEDIA_TYPE_VIDEO, -1, -1, nullptr, 0);
    if (index < 0) {
        throw VideoError("video: " + pathUtf8 + " has no video stream");
    }
    const AVCodecParameters* par = input.format->streams[index]->codecpar;
    const AVCodec* codec = avcodec_find_decoder(par->codec_id);
    if (codec == nullptr) {
        throw VideoError("video: no decoder for " + nameOr(avcodec_get_name(par->codec_id)) + " in this build");
    }
    DecoderCloser dec;
    dec.decoder = avcodec_alloc_context3(codec);
    dec.frame = av_frame_alloc();
    input.packet = av_packet_alloc();
    if (dec.decoder == nullptr || dec.frame == nullptr || input.packet == nullptr) {
        throw VideoError("video: out of memory");
    }
    avcodec_parameters_to_context(dec.decoder, par);
    dec.decoder->thread_count = 1;
    err = avcodec_open2(dec.decoder, codec, nullptr);
    if (err < 0) {
        throw VideoError("video: cannot open the decoder: " + detail::ffmpegError(err));
    }
    std::vector<std::vector<std::uint8_t>> frames;
    const auto receive = [&]() {
        for (;;) {
            const int r = avcodec_receive_frame(dec.decoder, dec.frame);
            if (r == AVERROR(EAGAIN) || r == AVERROR_EOF) return;
            if (r < 0) throw VideoError("video: cannot decode " + pathUtf8 + ": " + detail::ffmpegError(r));
            const int w = dec.frame->width, h = dec.frame->height;
            dec.scaler = sws_getCachedContext(dec.scaler, w, h, static_cast<AVPixelFormat>(dec.frame->format), w, h,
                                              AV_PIX_FMT_RGBA, SWS_BICUBIC | SWS_ACCURATE_RND | SWS_FULL_CHR_H_INT,
                                              nullptr, nullptr, nullptr);
            if (dec.scaler == nullptr) throw VideoError("video: cannot convert the decoded frames to RGBA");
            const int space = dec.frame->colorspace == AVCOL_SPC_BT709 ? SWS_CS_ITU709 : SWS_CS_ITU601;
            const int fullRange = dec.frame->color_range == AVCOL_RANGE_JPEG ? 1 : 0;
            sws_setColorspaceDetails(dec.scaler, sws_getCoefficients(space), fullRange, sws_getCoefficients(space), 1,
                                     0, 1 << 16, 1 << 16);
            std::vector<std::uint8_t> rgba(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 4);
            std::uint8_t* dst[4] = {rgba.data(), nullptr, nullptr, nullptr};
            const int dstStride[4] = {w * 4, 0, 0, 0};
            sws_scale(dec.scaler, dec.frame->data, dec.frame->linesize, 0, h, dst, dstStride);
            frames.push_back(std::move(rgba));
            av_frame_unref(dec.frame);
        }
    };
    while ((err = av_read_frame(input.format, input.packet)) >= 0) {
        if (input.packet->stream_index == index) {
            const int s = avcodec_send_packet(dec.decoder, input.packet);
            if (s < 0) {
                av_packet_unref(input.packet);
                throw VideoError("video: cannot decode " + pathUtf8 + ": " + detail::ffmpegError(s));
            }
            receive();
        }
        av_packet_unref(input.packet);
    }
    avcodec_send_packet(dec.decoder, nullptr);
    receive();
    return frames;
}

#else

VideoProbe probeVideo(const fs::path&) { throw VideoError("video encoding is not available in this build"); }
bool hasVideoDecoder(std::string_view) { return false; }
std::vector<std::vector<std::uint8_t>> decodeVideoFrames(const fs::path&) {
    throw VideoError("video encoding is not available in this build");
}

#endif

}  // namespace mdbridge::render::video
