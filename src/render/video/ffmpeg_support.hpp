#pragma once
// Internal helpers shared by encoder.cpp and probe.cpp (only compiled with the FFmpeg libraries).

#if MDBRIDGE_HAVE_FFMPEG

#include <string>

namespace mdbridge::render::video::detail {

// Once per process: FFmpeg's log goes to a small buffer instead of stderr (MDATK runs ffmpeg with
// `-loglevel error` and quotes its stderr in errors), and SVT-AV1's own log is limited to fatal
// errors unless SVT_LOG is set.
void initFfmpeg();
// The FFmpeg error messages logged since the last call, joined with "; " (empty when none).
std::string takeFfmpegLog();
// av_strerror, plus the logged messages: "Invalid argument (libsvtav1: ...)".
std::string ffmpegError(int code);

}  // namespace mdbridge::render::video::detail

#endif
