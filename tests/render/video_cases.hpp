#pragma once
// Helpers shared by the video frame loop tests (test_video_loop.cpp without Qt, test_video_qt.cpp
// with the QPainter compositor): the video runs of tests/data/render/cli/render_cli_cases.json
// (tests/parity/make_render_cli_cases.py: MDATK's `render` on the demo project, 24 frames at
// 320 x 180) as a project (the case's project file with its `-s` overrides) on the demo dataset
// (tests/data/render/inputs/demo_dry, the same analyses of the same trajectory), and the video
// MDATK wrote (ffprobe's view of it).

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "commands.hpp"
#include "figure_records.hpp"
#include "mdbridge/core/project_io.hpp"
#include "video/probe.hpp"

namespace mdbridge::test {

inline const core::Value& renderCliCases() {
    static const core::Value doc = readJson(kRenderData / "cli" / "render_cli_cases.json");
    return doc;
}

inline const core::Value& renderCliCase(std::string_view name) {
    for (const core::Value& c : *renderCliCases().find("cases")->get<core::Array>()) {
        if (c.find("name")->stringOr() == name) return c;
    }
    throw std::runtime_error("no render CLI case " + std::string(name));
}

// The case's project: its project file read back, then its `-s key=value` overrides.
inline core::Project renderCliProject(const core::Value& record) {
    TempDir dir;
    writeFile(dir / "demo.mdatk.toml", renderCliCases().find("project_file")->stringOr());
    core::Project project = core::loadProject(dir / "demo.mdatk.toml");
    const core::Array& argv = *record.find("argv")->get<core::Array>();
    for (std::size_t k = 0; k + 1 < argv.size(); ++k) {
        if (argv[k].stringOr() == "-s") app::applyOverride(project, argv[k + 1].stringOr());
    }
    return project;
}

inline const core::Dataset& demoDataset() {
    static const core::Dataset ds = loadRenderInput("demo_dry");
    return ds;
}

// What MDATK's video looked like (ffprobe), codec name excluded (RN-07: MDATK's is h264 or mpeg4,
// MDBridge's the bundled AV1 or mpeg4).
struct ExpectedVideo {
    std::string file;  // "results/MDATK_video_{stamp}.mp4"
    int width = 0, height = 0;
    std::string pixelFormat;
    int rateNum = 0, rateDen = 1;
    std::int64_t frames = 0;
    std::int64_t streams = 0, audioStreams = 0;
    std::string format;
    bool faststart = false;
};

inline ExpectedVideo expectedVideo(const core::Value& record) {
    for (const auto& [name, entry] : *record.find("files")->get<core::Table>()) {
        const core::Value* v = entry.find("video");
        if (!v) continue;
        const core::Value& s = (*v->find("streams")->get<core::Array>())[0];
        ExpectedVideo e;
        e.file = name;
        e.width = static_cast<int>(*s.find("width")->get<std::int64_t>());
        e.height = static_cast<int>(*s.find("height")->get<std::int64_t>());
        e.pixelFormat = s.find("pix_fmt")->stringOr();
        const std::string rate = s.find("r_frame_rate")->stringOr();
        e.rateNum = std::stoi(rate.substr(0, rate.find('/')));
        e.rateDen = std::stoi(rate.substr(rate.find('/') + 1));
        e.frames = std::stoll(s.find("nb_read_frames")->stringOr());
        e.streams = *v->find("format")->find("nb_streams")->get<std::int64_t>();
        e.format = v->find("format")->find("format_name")->stringOr();
        e.audioStreams = *v->find("audio_streams")->get<std::int64_t>();
        e.faststart = *v->find("faststart")->get<bool>();
        return e;
    }
    throw std::runtime_error("the case wrote no video");
}

// The case's stderr lines that start with `prefix` ("Rendering ", "Wrote "), with {work} replaced.
inline std::vector<std::string> caseLines(const core::Value& record, const std::string& prefix, const std::string& work) {
    std::vector<std::string> out;
    const std::string err = record.find("stderr")->stringOr();
    std::size_t at = 0;
    while (at < err.size()) {
        const std::size_t end = std::min(err.find('\n', at), err.size());
        std::string line = err.substr(at, end - at);
        at = end + 1;
        if (line.rfind(prefix, 0) != 0) continue;
        for (std::size_t p = line.find("{work}"); p != std::string::npos; p = line.find("{work}", p + work.size())) {
            line.replace(p, 6, work);
        }
        out.push_back(line);
    }
    return out;
}

// "{stamp}" -> a fixed time stamp in a case's file name
inline std::string stamped(std::string name) {
    const std::size_t p = name.find("{stamp}");
    if (p != std::string::npos) name.replace(p, 7, "2026-11-25_120000");
    return name;
}

// The written video against MDATK's: frame count, size, frame rate, pixel format, streams, no audio,
// container, faststart (not the codec name, RN-07).
inline void expectVideoLike(const render::video::VideoProbe& got, const ExpectedVideo& want) {
    EXPECT_EQ(got.packets, want.frames);
    EXPECT_EQ(got.nbFrames, want.frames);
    EXPECT_EQ(got.width, want.width);
    EXPECT_EQ(got.height, want.height);
    EXPECT_EQ(got.pixelFormat, want.pixelFormat);
    EXPECT_EQ(got.frameRateNum, want.rateNum);
    EXPECT_EQ(got.frameRateDen, want.rateDen);
    EXPECT_EQ(got.videoStreams + got.audioStreams + got.otherStreams, want.streams);
    EXPECT_EQ(got.audioStreams, want.audioStreams);
    EXPECT_EQ(got.format, want.format);
    EXPECT_EQ(got.moovBeforeMdat, want.faststart);
    // frame k at k / fps: the duration is frames / fps
    EXPECT_NEAR(got.durationSeconds, static_cast<double>(want.frames) * want.rateDen / want.rateNum, 1e-9);
}

}  // namespace mdbridge::test
