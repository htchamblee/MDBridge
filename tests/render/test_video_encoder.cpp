// The MP4 encoder (src/render/video): codec and quality mapping (RN-07), the LGPL-only check of
// the FFmpeg libraries in use, and encoded files checked by demuxing them (frame count, size,
// frame rate, time base, duration, codec, yuv420p, faststart, no audio), as the video parity tests
// compare them (docs/specs/m4-plan.md section 3). The encoding tests skip in builds without FFmpeg;
// CI builds with -DMDBRIDGE_FFMPEG=ON, so they run there.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include "mdbridge/core/paths.hpp"
#include "video/encoder.hpp"
#include "video/probe.hpp"

namespace fs = std::filesystem;
using namespace mdbridge::render::video;

namespace {

bool hasEncoder(const std::set<std::string>& encoders, std::string_view name) {
    return encoders.count(std::string(name)) != 0;
}

std::function<bool(std::string_view)> encoders(std::set<std::string> names) {
    return [names = std::move(names)](std::string_view name) { return hasEncoder(names, name); };
}

// A frame of a 24-frame test animation: colour gradients and a square moving across.
std::vector<std::uint8_t> testFrame(int width, int height, int index) {
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4);
    const int x0 = index * (width - 40) / 23;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            std::uint8_t* p = &rgba[(static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                                     static_cast<std::size_t>(x)) * 4];
            const bool square = x >= x0 && x < x0 + 40 && y >= height / 2 - 20 && y < height / 2 + 20;
            p[0] = square ? 255 : static_cast<std::uint8_t>(x * 255 / (width - 1));
            p[1] = square ? 255 : static_cast<std::uint8_t>(y * 255 / (height - 1));
            p[2] = square ? 255 : static_cast<std::uint8_t>(128 + index * 4);
            p[3] = 255;
        }
    }
    return rgba;
}

class VideoFiles : public ::testing::Test {
protected:
    void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        dir = fs::temp_directory_path() / ("mdbridge_video_" + std::string(info->test_suite_name()) + "_" + info->name());
        fs::remove_all(dir);
        fs::create_directories(dir);
    }
    void TearDown() override {
        std::error_code ec;
        fs::remove_all(dir, ec);
    }
    // Encodes `frames` test frames; returns the warnings.
    std::vector<std::string> encode(const fs::path& path, EncoderOptions options, int frames) {
        std::vector<std::string> warnings;
        VideoEncoder encoder(path, options, [&](const std::string& w) { warnings.push_back(w); });
        EXPECT_TRUE(fs::exists(encoder.temporaryPath()));
        EXPECT_EQ(encoder.temporaryPath().parent_path(), fs::absolute(path).parent_path());
        for (int k = 0; k < frames; ++k) {
            const std::vector<std::uint8_t> rgba = testFrame(options.width, options.height, k);
            encoder.writeFrame(rgba.data());
        }
        EXPECT_EQ(encoder.framesWritten(), frames);
        EXPECT_EQ(encoder.finish(), path);
        EXPECT_FALSE(fs::exists(encoder.temporaryPath()));
        return warnings;
    }
    // Nothing but `expected` in the folder (no temporary files left behind).
    void expectOnly(const std::set<std::string>& expected) {
        std::set<std::string> names;
        for (const auto& entry : fs::directory_iterator(dir)) {
            names.insert(entry.path().filename().string());
        }
        EXPECT_EQ(names, expected);
    }
    fs::path dir;
};

#define SKIP_WITHOUT_FFMPEG()                                                        \
    if (!videoEncodingAvailable()) {                                                 \
        GTEST_SKIP() << "built without the FFmpeg libraries (packaging/ffmpeg)"; \
    }

}  // namespace

// ---- codec choice and quality (no FFmpeg needed) ----

TEST(VideoCodecChoice, H264AndH265RequestsUseTheBundledAv1Encoder) {
    const auto both = encoders({"libsvtav1", "mpeg4"});
    CodecChoice c = chooseCodec("libx264", 20, both);
    EXPECT_EQ(c.requested, "libx264");
    EXPECT_EQ(c.encoder, "libsvtav1");
    EXPECT_EQ(c.message, "video: codec libx264 is not bundled (GPL); using libsvtav1 (AV1)");
    EXPECT_EQ(c.qualityOption, "crf");
    EXPECT_EQ(c.quality, 25);
    EXPECT_EQ(chooseCodec("libx265", 28, both).message, "video: codec libx265 is not bundled (GPL); using libsvtav1 (AV1)");
    EXPECT_EQ(chooseCodec("libx265", 28, both).quality, 35);
    EXPECT_EQ(chooseCodec("h264", 20, both).message, "video: codec h264 is not bundled; using libsvtav1 (AV1)");
    EXPECT_EQ(chooseCodec("hevc", 20, both).encoder, "libsvtav1");
    // the project files' default when the codec is empty
    EXPECT_EQ(chooseCodec("", 20, both).requested, "libx264");
    EXPECT_EQ(chooseCodec("", 20, both).encoder, "libsvtav1");
}

TEST(VideoCodecChoice, BundledEncodersAreUsedAsNamedWithoutMessage) {
    const auto both = encoders({"libsvtav1", "mpeg4"});
    const CodecChoice av1 = chooseCodec("libsvtav1", 30, both);
    EXPECT_EQ(av1.encoder, "libsvtav1");
    EXPECT_EQ(av1.message, "");
    EXPECT_EQ(av1.quality, 37);
    EXPECT_EQ(chooseCodec("av1", 20, both).encoder, "libsvtav1");
    EXPECT_EQ(chooseCodec("av1", 20, both).message, "");
    const CodecChoice mpeg4 = chooseCodec("mpeg4", 20, both);
    EXPECT_EQ(mpeg4.encoder, "mpeg4");
    EXPECT_EQ(mpeg4.message, "");
    EXPECT_EQ(mpeg4.qualityOption, "qscale");
    EXPECT_EQ(mpeg4.quality, 3);  // MDATK's `-q:v 3` at the default crf
}

TEST(VideoCodecChoice, OtherNamesPassThroughOrFallBack) {
    // an encoder the libraries have is used as named, with its own default quality (as MDATK)
    const CodecChoice other = chooseCodec("h264_mf", 20, encoders({"libsvtav1", "mpeg4", "h264_mf"}));
    EXPECT_EQ(other.encoder, "h264_mf");
    EXPECT_EQ(other.message, "");
    EXPECT_EQ(other.qualityOption, "");
    // one they lack falls back to the first of libsvtav1 / mpeg4 with a warning
    const auto both = encoders({"libsvtav1", "mpeg4"});
    EXPECT_EQ(chooseCodec("libopenh264", 20, both).message, "video: encoder libopenh264 is not available; using libsvtav1");
    EXPECT_EQ(chooseCodec("libaom-av1", 20, both).encoder, "libsvtav1");
    const auto onlyMpeg4 = encoders({"mpeg4"});
    const CodecChoice fallback = chooseCodec("libx264", 20, onlyMpeg4);
    EXPECT_EQ(fallback.encoder, "mpeg4");
    EXPECT_EQ(fallback.message, "video: encoder libx264 is not available; using mpeg4");
    EXPECT_EQ(fallback.quality, 3);
    EXPECT_EQ(chooseCodec("libsvtav1", 20, onlyMpeg4).message, "video: encoder libsvtav1 is not available; using mpeg4");
    try {
        chooseCodec("libx264", 20, encoders({}));
        FAIL() << "no exception";
    } catch (const VideoError& e) {
        EXPECT_STREQ(e.what(), "no video encoder is available (tried libx264, libsvtav1 and mpeg4)");
        EXPECT_EQ(e.pythonName(), "VideoError");
    }
}

TEST(VideoCodecChoice, CrfTables) {
    // x264 / x265 crf -> SVT-AV1 crf: round(crf * 63 / 51), clamped to 0..51 first
    const std::vector<std::pair<int, int>> svt = {{-5, 0}, {0, 0}, {1, 1}, {10, 12}, {18, 22}, {20, 25},
                                                   {23, 28}, {28, 35}, {35, 43}, {51, 63}, {70, 63}};
    for (const auto& [crf, expected] : svt) {
        EXPECT_EQ(svtAv1Crf(crf), expected) << crf;
    }
    // -> mpeg4 qscale: min(31, max(2, round(3 * 2 ** ((crf - 20) / 6)))) of MDATK 0.1.2 (Python's
    // values for crf 0 ... 51; 20 -> MDATK 0.0.4.1's q 3)
    const std::vector<int> q = {2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                                2,  3,  3,  3,  4,  4,  5,  5,  6,  7,  8,  8,  10, 11, 12, 13, 15, 17,
                                19, 21, 24, 27, 30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31};
    for (int crf = 0; crf <= 51; ++crf) {
        EXPECT_EQ(mpeg4Qscale(crf), q[static_cast<std::size_t>(crf)]) << crf;
    }
    EXPECT_EQ(mpeg4Qscale(-1), 2);
    EXPECT_EQ(mpeg4Qscale(99), 31);
    for (int crf = 1; crf <= 51; ++crf) {  // monotonic
        EXPECT_GE(svtAv1Crf(crf), svtAv1Crf(crf - 1));
        EXPECT_GE(mpeg4Qscale(crf), mpeg4Qscale(crf - 1));
    }
}

TEST(VideoCodecChoice, TemporaryNameAsMdatk) {
    const fs::path tmp = temporaryVideoPath(fs::path("out") / "MDATK_video.mp4");
    EXPECT_EQ(tmp.parent_path(), fs::absolute(fs::path("out")));
    const std::string name = tmp.filename().string();
    EXPECT_EQ(name.rfind(".MDATK_tmp_", 0), 0u) << name;
    EXPECT_EQ(name.substr(name.size() - std::string("_MDATK_video.mp4").size()), "_MDATK_video.mp4") << name;
}

// ---- the FFmpeg libraries ----

TEST(VideoFfmpeg, LgplOnlyBuild) {
    if (!videoEncodingAvailable()) {
        EXPECT_THROW(ffmpegBuildInfo(), VideoError);
        EXPECT_EQ(ffmpegVersionLine(), "not available");
        GTEST_SKIP() << "built without the FFmpeg libraries";
    }
    const FfmpegBuildInfo info = ffmpegBuildInfo();
    EXPECT_TRUE(info.lgplOnly) << info.license << ": " << info.configuration;
    EXPECT_EQ(info.license, "LGPL version 2.1 or later");
    EXPECT_EQ(info.configuration.find("--enable-gpl"), std::string::npos) << info.configuration;
    EXPECT_EQ(info.configuration.find("--enable-nonfree"), std::string::npos) << info.configuration;
    EXPECT_NE(info.configuration.find("--enable-libsvtav1"), std::string::npos) << info.configuration;
    EXPECT_TRUE(info.hasSvtAv1);
    EXPECT_TRUE(info.hasMpeg4);
    EXPECT_FALSE(info.version.empty());
    const std::string line = ffmpegVersionLine();
    EXPECT_EQ(line.rfind("FFmpeg " + info.version + " libraries (LGPL version 2.1 or later): libavcodec ", 0), 0u)
        << line;
}

TEST(VideoFfmpeg, NotAvailableMessage) {
    if (videoEncodingAvailable()) {
        GTEST_SKIP() << "built with the FFmpeg libraries";
    }
    try {
        VideoEncoder encoder("x.mp4", EncoderOptions{320, 180, 30.0, "libx264", 20, 0});
        FAIL() << "no exception";
    } catch (const VideoError& e) {
        EXPECT_STREQ(e.what(), "video encoding is not available in this build");
    }
}

// ---- encoded files ----

TEST_F(VideoFiles, Av1Default24Frames320x180) {
    SKIP_WITHOUT_FFMPEG();
    const fs::path path = dir / "MDATK_video.mp4";
    const auto warnings = encode(path, EncoderOptions{320, 180, 30.0, "libx264", 20, 0}, 24);
    ASSERT_EQ(warnings, std::vector<std::string>{"video: codec libx264 is not bundled (GPL); using libsvtav1 (AV1)"});
    expectOnly({"MDATK_video.mp4"});

    const VideoProbe p = probeVideo(path);
    EXPECT_EQ(p.format, "mov,mp4,m4a,3gp,3g2,mj2");
    EXPECT_EQ(p.videoStreams, 1);
    EXPECT_EQ(p.audioStreams, 0);
    EXPECT_EQ(p.otherStreams, 0);
    EXPECT_EQ(p.codec, "av1");
    EXPECT_EQ(p.pixelFormat, "yuv420p");
    EXPECT_EQ(p.width, 320);
    EXPECT_EQ(p.height, 180);
    EXPECT_EQ(p.frameRateNum, 30);
    EXPECT_EQ(p.frameRateDen, 1);
    // what ffmpeg -r 30 writes into MP4: time base 1/15360 (1/30 doubled up to >= 10000)
    EXPECT_EQ(p.timeBaseNum, 1);
    EXPECT_EQ(p.timeBaseDen, 15360);
    EXPECT_EQ(p.durationTs, 24 * 512);
    EXPECT_DOUBLE_EQ(p.durationSeconds, 0.8);
    EXPECT_EQ(p.nbFrames, 24);
    EXPECT_EQ(p.packets, 24);
    EXPECT_GE(p.keyframes, 1);
    EXPECT_EQ(p.firstPts, 0);
    EXPECT_TRUE(p.moovBeforeMdat);
    EXPECT_EQ(p.colorSpace, "bt709");
    EXPECT_EQ(p.colorRange, "tv");
    EXPECT_EQ(p.colorPrimaries, "bt709");
    EXPECT_EQ(p.colorTransfer, "bt709");
}

TEST_F(VideoFiles, Mpeg4Fallback24Frames320x180) {
    SKIP_WITHOUT_FFMPEG();
    const fs::path path = dir / "clip.mp4";
    const auto warnings = encode(path, EncoderOptions{320, 180, 24.0, "mpeg4", 20, 0}, 24);
    EXPECT_TRUE(warnings.empty());
    const VideoProbe p = probeVideo(path);
    EXPECT_EQ(p.codec, "mpeg4");
    EXPECT_EQ(p.pixelFormat, "yuv420p");
    EXPECT_EQ(p.width, 320);
    EXPECT_EQ(p.height, 180);
    EXPECT_EQ(p.frameRateNum, 24);
    EXPECT_EQ(p.frameRateDen, 1);
    EXPECT_EQ(p.timeBaseNum, 1);
    EXPECT_EQ(p.timeBaseDen, 12288);
    EXPECT_EQ(p.durationTs, 24 * 512);
    EXPECT_DOUBLE_EQ(p.durationSeconds, 1.0);
    EXPECT_EQ(p.packets, 24);
    EXPECT_EQ(p.audioStreams, 0);
    EXPECT_TRUE(p.moovBeforeMdat);
}

// The RGBA -> yuv420p conversion runs on several threads (EncoderOptions::conversionThreads): every
// thread computes whole rows from the complete frame, so the video is the same bytes for any number
// of threads; rows passed with padding (a stride) give the same video as packed rows.
TEST_F(VideoFiles, ConversionThreadsGiveTheSameVideo) {
    SKIP_WITHOUT_FFMPEG();
    constexpr int w = 640, h = 360, frames = 5;
    std::vector<std::vector<std::uint8_t>> input;
    std::uint32_t state = 12345u;
    for (int k = 0; k < frames; ++k) {
        std::vector<std::uint8_t> rgba = testFrame(w, h, k * 5);
        for (std::size_t i = 0; i < rgba.size(); i += 4) {  // noise: chroma filtering sees every row
            for (int c = 0; c < 3; ++c) {
                state = state * 1664525u + 1013904223u;
                rgba[i + static_cast<std::size_t>(c)] ^= static_cast<std::uint8_t>((state >> 24) & 0x1F);
            }
        }
        input.push_back(std::move(rgba));
    }
    const auto encodeWith = [&](int conversionThreads, std::size_t padding) {
        EncoderOptions options{w, h, 24.0, "mpeg4", 20, 1};
        options.conversionThreads = conversionThreads;
        const fs::path path = dir / ("t" + std::to_string(conversionThreads) + "_" + std::to_string(padding) + ".mp4");
        VideoEncoder encoder(path, options);
        const std::size_t stride = static_cast<std::size_t>(w) * 4 + padding;
        for (const std::vector<std::uint8_t>& rgba : input) {
            std::vector<std::uint8_t> padded(stride * h, 0x5A);
            for (int y = 0; y < h; ++y) {
                std::copy_n(rgba.begin() + static_cast<std::ptrdiff_t>(y) * w * 4, w * 4,
                            padded.begin() + static_cast<std::ptrdiff_t>(stride * static_cast<std::size_t>(y)));
            }
            encoder.writeFrame(padded.data(), stride);
        }
        encoder.finish();
        std::ifstream in(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
    };
    const std::string one = encodeWith(1, 0);
    ASSERT_GT(one.size(), 1000u);
    EXPECT_TRUE(one == encodeWith(3, 0));
    EXPECT_TRUE(one == encodeWith(8, 0));
    EXPECT_TRUE(one == encodeWith(4, 64));
    if (hasVideoDecoder("mpeg4")) {
        EXPECT_EQ(decodeVideoFrames(dir / "t1_0.mp4").size(), static_cast<std::size_t>(frames));
    }
}

TEST_F(VideoFiles, FractionalFrameRateAndReplacingAFile) {
    SKIP_WITHOUT_FFMPEG();
    const fs::path path = dir / "video.mp4";
    {
        std::ofstream(path) << "old";
    }
    encode(path, EncoderOptions{64, 48, 29.97, "libsvtav1", 20, 0}, 6);
    expectOnly({"video.mp4"});
    const VideoProbe p = probeVideo(path);
    EXPECT_EQ(p.codec, "av1");
    EXPECT_EQ(p.frameRateNum, 2997);  // av_d2q(29.97, 1001000), as `ffmpeg -r 29.97`
    EXPECT_EQ(p.frameRateDen, 100);
    EXPECT_EQ(p.packets, 6);
}

TEST_F(VideoFiles, NonAsciiFolder) {
    SKIP_WITHOUT_FFMPEG();  // FFmpeg's file protocol takes UTF-8 names (wide-character API on Windows)
    const fs::path folder = dir / mdbridge::core::pathFromUtf8("vid\xc3\xa9os \xce\xb1\xce\xb2");
    fs::create_directories(folder);
    const fs::path path = folder / "MDATK_video.mp4";
    encode(path, EncoderOptions{64, 36, 30.0, "libx264", 20, 0}, 4);
    EXPECT_TRUE(fs::exists(path));
    EXPECT_EQ(probeVideo(path).packets, 4);
    EXPECT_TRUE(probeVideo(path).moovBeforeMdat);
}

TEST_F(VideoFiles, SizesSvtAv1CannotEncodeFallBackToMpeg4) {
    SKIP_WITHOUT_FFMPEG();
    const fs::path path = dir / "tiny.mp4";
    const auto warnings = encode(path, EncoderOptions{2, 2, 30.0, "libx264", 20, 0}, 3);
    ASSERT_EQ(warnings.size(), 2u);
    EXPECT_EQ(warnings[0], "video: codec libx264 is not bundled (GPL); using libsvtav1 (AV1)");
    EXPECT_EQ(warnings[1].rfind("video: libsvtav1 cannot encode 2x2 (", 0), 0u) << warnings[1];
    EXPECT_EQ(warnings[1].find(" @ 0x"), std::string::npos) << warnings[1];
    EXPECT_EQ(warnings[1].substr(warnings[1].size() - 14), "); using mpeg4") << warnings[1];
    const VideoProbe p = probeVideo(path);
    EXPECT_EQ(p.codec, "mpeg4");
    EXPECT_EQ(p.width, 2);
    EXPECT_EQ(p.packets, 3);
}

TEST_F(VideoFiles, AbortAndDestructorRemoveTheTemporaryFile) {
    SKIP_WITHOUT_FFMPEG();
    const fs::path path = dir / "MDATK_video.mp4";
    const EncoderOptions options{320, 180, 30.0, "libx264", 20, 0};
    {
        VideoEncoder encoder(path, options);
        const std::vector<std::uint8_t> rgba = testFrame(320, 180, 0);
        encoder.writeFrame(rgba.data());
        EXPECT_TRUE(fs::exists(encoder.temporaryPath()));
        encoder.abort();
        EXPECT_FALSE(fs::exists(encoder.temporaryPath()));
        EXPECT_THROW(encoder.writeFrame(rgba.data()), VideoError);
        EXPECT_THROW(encoder.finish(), VideoError);
    }
    {
        VideoEncoder encoder(path, options);  // e.g. cancelled: destroyed without finish()
        const std::vector<std::uint8_t> rgba = testFrame(320, 180, 0);
        encoder.writeFrame(rgba.data());
    }
    expectOnly({});
}

TEST_F(VideoFiles, ClearErrors) {
    SKIP_WITHOUT_FFMPEG();
    try {
        VideoEncoder encoder(dir / "a.mp4", EncoderOptions{321, 180, 30.0, "libx264", 20, 0});
        FAIL() << "no exception";
    } catch (const VideoError& e) {
        EXPECT_STREQ(e.what(), "video size 321x180 is not even (yuv420p needs even sizes of at least 2x2)");
    }
    try {
        VideoEncoder encoder(dir / "a.mp4", EncoderOptions{320, 180, 0.0, "libx264", 20, 0});
        FAIL() << "no exception";
    } catch (const VideoError& e) {
        EXPECT_EQ(std::string(e.what()).rfind("video frame rate 0", 0), 0u) << e.what();
    }
    try {
        VideoEncoder encoder(dir / "missing" / "a.mp4", EncoderOptions{320, 180, 30.0, "mpeg4", 20, 0});
        FAIL() << "no exception";
    } catch (const VideoError& e) {
        EXPECT_EQ(std::string(e.what()).rfind("video: cannot write ", 0), 0u) << e.what();
    }
    {
        VideoEncoder encoder(dir / "a.mp4", EncoderOptions{320, 180, 30.0, "mpeg4", 20, 0});
        EXPECT_THROW(encoder.finish(), VideoError);  // no frames
    }
    expectOnly({});
}
