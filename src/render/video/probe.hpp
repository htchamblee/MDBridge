#pragma once
// Self-check of a written video by demuxing it (no decoding): what ffprobe would report about the
// first video stream, as the parity tests compare it (docs/specs/m4-plan.md section 3, "Video":
// frame count, size, r_frame_rate, time base and duration, codec, yuv420p, moov before mdat, no
// audio). Needs the FFmpeg libraries (VideoError otherwise); the box order is read directly.

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::render::video {

struct VideoProbe {
    std::string format;     // demuxer, e.g. "mov,mp4,m4a,3gp,3g2,mj2"
    int videoStreams = 0;
    int audioStreams = 0;
    int otherStreams = 0;
    // the first video stream
    std::string codec;      // avcodec_get_name: "av1", "mpeg4", ...
    std::string pixelFormat;  // "yuv420p" (from the av1 / mpeg4 parser or decoder)
    int width = 0;
    int height = 0;
    int frameRateNum = 0;   // r_frame_rate
    int frameRateDen = 1;
    int timeBaseNum = 0;    // the stream's time base
    int timeBaseDen = 1;
    std::int64_t durationTs = 0;  // stream duration in time-base units
    double durationSeconds = 0.0;
    std::int64_t nbFrames = 0;    // the container's sample count
    std::int64_t packets = 0;     // packets demuxed (the real frame count)
    std::int64_t keyframes = 0;
    std::int64_t firstPts = 0;
    std::string colorSpace;       // av_color_space_name, e.g. "bt709" ("unknown" when untagged)
    std::string colorRange;       // "tv" / "pc" / "unknown"
    std::string colorPrimaries;
    std::string colorTransfer;
    // MP4 box order (top level): faststart puts moov before mdat
    bool moovBeforeMdat = false;
};

VideoProbe probeVideo(const std::filesystem::path& path);

// The FFmpeg libraries have a decoder for this codec ("mpeg4": yes; "av1": not in packaging/ffmpeg's
// build, which only encodes AV1). False without FFmpeg.
bool hasVideoDecoder(std::string_view codec);

// The decoded frames of the first video stream in order, as straight RGBA (width x height x 4, rows
// top first, alpha 255), converted from yuv420p with the matrix and range the stream is tagged with
// (MDBridge's videos: BT.709, limited range) - the inverse of the encoder's conversion. For the
// tests' PSNR check of encoded frames against the composed ones. Throws VideoError when there is
// no decoder for the stream's codec (hasVideoDecoder) or the file cannot be decoded.
std::vector<std::vector<std::uint8_t>> decodeVideoFrames(const std::filesystem::path& path);

}  // namespace mdbridge::render::video
