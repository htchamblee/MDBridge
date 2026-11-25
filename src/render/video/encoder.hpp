#pragma once
// MP4 video encoding with the bundled LGPL FFmpeg libraries (MDATK: render/video.py, FFmpegWriter).
//
// MDATK pipes raw RGBA frames into an `ffmpeg` executable; MDBridge links FFmpeg's libraries
// (libavcodec, libavformat, libavutil, libswscale; packaging/ffmpeg) and does the same in process:
// RGBA -> yuv420p (swscale), the chosen encoder, an MP4 file with the `moov` box in front
// (`-movflags +faststart`), no audio, frame k at time k / fps. The file is written under a temporary
// name in the target folder (`.MDATK_tmp_<pid>_<name>`, as MDATK) and renamed when finished;
// failures and abort() remove it.
//
// Codecs (RN-07, docs/DIFFERENCES.md "Video encoding"): `libx264` (the project files' default),
// `libx265`, `h264`, `hevc` and `av1` map to the bundled AV1 encoder (SVT-AV1); `libsvtav1` and
// `mpeg4` are used as named; any other name is used if the FFmpeg libraries have such an encoder,
// else the first available of libsvtav1 / mpeg4 with a warning. `crf` is on the x264 / x265 scale
// (0-51, default 20) and is mapped to the encoder's own quality scale (svtAv1Crf, mpeg4Qscale).
//
// Without FFmpeg (a source build that did not find packaging/ffmpeg's libraries), every entry point
// that needs the libraries throws VideoError("video encoding is not available in this build").

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <string_view>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::render::video {

// MDATK's VideoError (a RuntimeError): step-log summaries read "VideoError: <message>".
class VideoError : public core::NamedError {
public:
    explicit VideoError(const std::string& message) : core::NamedError("VideoError", message) {}
};

// True when MDBridge was built with the FFmpeg libraries.
bool videoEncodingAvailable();

// What the FFmpeg libraries in use report about themselves (run time, not build time).
struct FfmpegBuildInfo {
    std::string version;        // FFmpeg release, e.g. "9.0.2" (av_version_info())
    std::string configuration;  // avcodec_configuration(): the configure flags
    std::string license;        // avcodec_license(), e.g. "LGPL version 2.1 or later"
    std::string libraries;      // "libavcodec 62.11.100, libavformat ..., libavutil ..., libswscale ..."
    bool lgplOnly = false;      // no --enable-gpl / --enable-nonfree / --enable-version3, LGPL 2.1+ licence
    bool hasSvtAv1 = false;     // the libsvtav1 encoder is present
    bool hasMpeg4 = false;      // the mpeg4 encoder is present
};
FfmpegBuildInfo ffmpegBuildInfo();  // throws VideoError without FFmpeg

// One line for the params log's `environment` (MDATK records the `ffmpeg -version` line), e.g.
// "FFmpeg 9.0.2 libraries (LGPL version 2.1 or later): libavcodec 62.11.100, ...", or
// "not available" without FFmpeg.
std::string ffmpegVersionLine();

// ---- codec choice (no FFmpeg needed; `available` says which encoders the libraries have) ----

inline constexpr std::string_view kAv1Encoder = "libsvtav1";
inline constexpr std::string_view kFallbackEncoder = "mpeg4";

// x264 / x265 CRF (clamped to 0..51) -> SVT-AV1 CRF (0..63): round(crf * 63 / 51); 20 -> 25.
int svtAv1Crf(int crf);
// x264 / x265 CRF (clamped to 0..51) -> mpeg4 qscale (1..31, `-q:v`): the default crf 20 gives
// MDATK's fixed `-q:v 3`; 1 + crf / 10 below 20, 3 + (crf - 20) * 28 / 31 above (rounded).
int mpeg4Qscale(int crf);

struct CodecChoice {
    std::string requested;  // the project's codec ("" -> "libx264")
    std::string encoder;    // the FFmpeg encoder to use
    // "video: codec libx264 is not bundled (GPL); using libsvtav1 (AV1)" and the like, or empty
    // when the encoder is the one requested
    std::string message;
    std::string qualityOption;  // "crf" (libsvtav1), "qscale" (mpeg4) or "" (encoder default)
    int quality = 0;            // its value
};
// Throws VideoError when no usable encoder is available.
CodecChoice chooseCodec(std::string_view requested, int crf,
                        const std::function<bool(std::string_view)>& available);

// ---- the encoder ----

struct EncoderOptions {
    int width = 0;   // even, >= 2 (yuv420p); MDATK's figure_size() makes them so
    int height = 0;
    double fps = 30.0;  // > 0; the frame rate as FFmpeg's `-r` reads it (av_d2q, 1001000)
    std::string codec = "libx264";
    int crf = 20;
    int threads = 0;  // encoder threads, 0 = the encoder's choice (all cores)
};

class VideoEncoder {
public:
    // Opens `.MDATK_tmp_<pid>_<name>` next to `path` (the folder must exist) and the encoder.
    // `warn` receives the codec message (CodecChoice::message) and other one-line warnings.
    VideoEncoder(const std::filesystem::path& path, const EncoderOptions& options,
                 std::function<void(const std::string&)> warn = {});
    ~VideoEncoder();  // abort() unless finish() succeeded
    VideoEncoder(const VideoEncoder&) = delete;
    VideoEncoder& operator=(const VideoEncoder&) = delete;

    // One frame: width x height straight (not premultiplied) RGBA, 8 bits per channel, rows
    // `strideBytes` apart (0 = width * 4); alpha is ignored, as FFmpeg's rgba -> yuv420p does.
    void writeFrame(const std::uint8_t* rgba, std::size_t strideBytes = 0);
    // Flushes the encoder, writes the trailer (faststart: moov moved in front) and renames the
    // temporary file to the target. Returns the target path.
    std::filesystem::path finish();
    // Stops and removes the temporary file; never throws.
    void abort() noexcept;

    int framesWritten() const;
    const CodecChoice& codec() const;
    const std::filesystem::path& temporaryPath() const;
    const std::filesystem::path& path() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m;
};

// The temporary name MDATK uses: `<folder>/.MDATK_tmp_<pid>_<file name>`.
std::filesystem::path temporaryVideoPath(const std::filesystem::path& path);

}  // namespace mdbridge::render::video
