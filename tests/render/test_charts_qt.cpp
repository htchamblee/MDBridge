// The standalone charts written by the Qt backend (render::writeFigure, figure_qt.hpp, and render
// core's PostScript / TIFF writers): for the chart cases of every format (data/render/charts,
// tests/parity/make_chart_cases.py) the files MDATK wrote are written again, each starting with the
// bytes MDATK's started with (PNG, JPEG / JFIF, TIFF little endian, SVG's XML declaration, PDF,
// PostScript) and of the same size (pixels: exact; SVG / PDF / EPS / PS pages: within a point, the
// figure being whole pixels). Contents are test_charts.cpp's; owner decision D10: no reference
// images. Runs on the offscreen platform (no display).

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <QColor>
#include <QGuiApplication>
#include <QImage>

#include "figure_records.hpp"
#include "mdbridge/render/figure_formats.hpp"
#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/panels/charts.hpp"
#include "mdbridge/render/text_qt.hpp"

using namespace mdbridge;
using namespace mdbridge::render;
using core::Array;
using core::Value;

namespace {

const std::filesystem::path kCharts = test::kRenderData / "charts";

std::string str(const Value& v) { return v.stringOr(); }

std::string bytesOf(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

std::string hex8(const std::string& data) {
    static const char digits[] = "0123456789abcdef";
    std::string out;
    for (std::size_t i = 0; i < 8 && i < data.size(); ++i) {
        const auto c = static_cast<unsigned char>(data[i]);
        out.push_back(digits[c >> 4]);
        out.push_back(digits[c & 15]);
    }
    return out;
}

std::uint32_t be32(const std::string& d, std::size_t at) {
    return (std::uint32_t{static_cast<unsigned char>(d[at])} << 24) | (std::uint32_t{static_cast<unsigned char>(d[at + 1])} << 16) |
           (std::uint32_t{static_cast<unsigned char>(d[at + 2])} << 8) | std::uint32_t{static_cast<unsigned char>(d[at + 3])};
}

// [width, height] in pixels (PNG, JPEG, TIFF) or points (SVG, PDF, EPS / PS) as the generator reads them
std::vector<double> sizeOf(const std::string& d) {
    std::smatch m;
    if (d.rfind("\x89PNG", 0) == 0) return {double(be32(d, 16)), double(be32(d, 20))};
    if (d.rfind("\xff\xd8", 0) == 0) {
        const QImage img = QImage::fromData(QByteArray(d.data(), static_cast<qsizetype>(d.size())), "JPEG");
        return {double(img.width()), double(img.height())};
    }
    if (d.rfind(std::string("II*\0", 4), 0) == 0) {
        auto u16 = [&](std::size_t at) { return unsigned(static_cast<unsigned char>(d[at])) | unsigned(static_cast<unsigned char>(d[at + 1])) << 8; };
        auto u32 = [&](std::size_t at) { return u16(at) | u16(at + 2) << 16; };
        const unsigned ifd = u32(4), n = u16(ifd);
        double w = 0, h = 0;
        for (unsigned k = 0; k < n; ++k) {
            const std::size_t e = ifd + 2 + 12 * k;
            if (u16(e) == 256) w = u32(e + 8);
            if (u16(e) == 257) h = u32(e + 8);
        }
        return {w, h};
    }
    if (std::regex_search(d, m, std::regex("<svg[^>]*width=\"([0-9.]+)mm\" height=\"([0-9.]+)mm\""))) {
        return {std::stod(m[1]) / 25.4 * 72.0, std::stod(m[2]) / 25.4 * 72.0};
    }
    if (std::regex_search(d, m, std::regex("/MediaBox \\[0 0 ([0-9.]+) ([0-9.]+)\\]"))) return {std::stod(m[1]), std::stod(m[2])};
    if (std::regex_search(d, m, std::regex("%%BoundingBox: 0 0 ([0-9]+) ([0-9]+)"))) return {std::stod(m[1]), std::stod(m[2])};
    return {};
}

class ChartFiles : public ::testing::TestWithParam<std::string> {};

TEST_P(ChartFiles, AsMdatkWroteThem) {
    const Value rec = test::readJson(kCharts / (GetParam() + ".json"));
    const core::Dataset ds = test::loadRenderInput(str(*rec.find("input")));
    const core::Project project = test::caseProject(rec);
    const std::filesystem::path tmp = std::filesystem::temp_directory_path() / ("mdbridge-charts-qt-" + GetParam());
    std::filesystem::remove_all(tmp);
    const charts::FigureWriter writer = [](const Figure& fig, const std::filesystem::path& path, const std::string& format) {
        writeFigure(fig, path, format);
    };
    const std::vector<std::string> paths =
        charts::renderCharts(ds, project, tmp / "out", core::Reporter(), str(*rec.find("stamp")), bundledFontSet(), writer);
    std::size_t k = 0;
    for (const Value& chart : *rec.find("charts")->get<Array>()) {
        const std::vector<double> face = test::numbersOf(*chart.path({"figure", "facecolor"}));
        for (const Value& file : *chart.find("files")->get<Array>()) {
            ASSERT_LT(k, paths.size());
            const std::string format = str(*file.find("format"));
            SCOPED_TRACE(paths[k]);
            const std::string data = bytesOf(paths[k++]);
            EXPECT_EQ(hex8(data), str(*file.find("magic"))) << format;
            const Value* px = file.find("size_px");
            const Value* pt = file.find("size_pt");
            const std::vector<double> size = sizeOf(data);
            ASSERT_EQ(size.size(), 2u) << format;
            if (px && !px->isNull()) {
                const std::vector<double> want = test::numbersOf(*px);
                EXPECT_EQ(size, want) << format;
            } else if (pt && !pt->isNull()) {
                std::vector<double> want = test::numbersOf(*pt);
                if (want.size() == 4) want = {want[2] - want[0], want[3] - want[1]};  // a MediaBox
                EXPECT_NEAR(size[0], want[0], 1.0) << format;
                EXPECT_NEAR(size[1], want[1], 1.0) << format;
            }
            if (format == "png") {
                // the figure's face colour in the corner
                const QImage img = QImage::fromData(QByteArray(data.data(), static_cast<qsizetype>(data.size())), "PNG");
                const QColor c = img.pixelColor(0, 0);
                EXPECT_NEAR(c.redF(), face[0], 1.5 / 255);
                EXPECT_NEAR(c.greenF(), face[1], 1.5 / 255);
                EXPECT_NEAR(c.blueF(), face[2], 1.5 / 255);
            }
        }
    }
    EXPECT_EQ(k, paths.size());
    if (!std::getenv("MDBRIDGE_KEEP_CHART_FILES")) std::filesystem::remove_all(tmp);  // kept: to look at them
}

INSTANTIATE_TEST_SUITE_P(Charts, ChartFiles,
                         ::testing::Values("formats", "format_jpg_eps", "format_tif_ps", "heatmaps", "dark",
                                           "dark_heatmaps", "style_default_smoothing", "stamped"),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

TEST(ChartPostScript, VectorsTextAndImage) {
    // a matrix chart as EPS: a valid DSC header, the page of the figure, text as outlines (no
    // fonts), the image's hex data, balanced graphics states
    const core::Dataset ds = test::loadRenderInput("demo_dry");
    core::Project p;
    p.output.charts.dpi = 40;
    const Figure fig = charts::matrixChart(ds, p, "contacts/per_residue_a", bundledFontSet());
    const std::string eps = renderPostScript(buildDisplayList(fig, bundledFontSet()), bundledFontSet(), true, "chart");
    EXPECT_EQ(eps.rfind("%!PS-Adobe-3.0 EPSF-3.0\n", 0), 0u);
    EXPECT_NE(eps.find("%%BoundingBox: 0 0 576 324\n"), std::string::npos);
    EXPECT_NE(eps.find("%%Title: chart\n"), std::string::npos);
    EXPECT_NE(eps.find("/ASCIIHexDecode filter >> image"), std::string::npos);
    EXPECT_EQ(eps.find("findfont"), std::string::npos);
    EXPECT_NE(eps.find(" c\n"), std::string::npos);  // glyph curves
    std::size_t saves = 0, restores = 0;
    for (std::size_t at = eps.find("gsave"); at != std::string::npos; at = eps.find("gsave", at + 1)) ++saves;
    for (std::size_t at = eps.find("grestore"); at != std::string::npos; at = eps.find("grestore", at + 1)) ++restores;
    EXPECT_EQ(saves, restores);
    EXPECT_EQ(eps.substr(eps.size() - 6), "%%EOF\n");
    const std::string ps = renderPostScript(buildDisplayList(fig, bundledFontSet()), bundledFontSet(), false);
    EXPECT_EQ(ps.rfind("%!PS-Adobe-3.0\n", 0), 0u);
    EXPECT_NE(ps.find("<< /PageSize [576 324] >> setpagedevice"), std::string::npos);
}

TEST(ChartTiff, BaselineRgba) {
    // read back by its tags (Qt's TIFF plugin is optional; pillow reads these files too)
    const std::vector<std::uint8_t> px{255, 0, 0, 255, 0, 255, 0, 128, 0, 0, 255, 0, 10, 20, 30, 40};
    const std::string tiff = encodeTiff(2, 2, px, 300.0);
    EXPECT_EQ(hex8(tiff), "49492a0008000000");  // as pillow writes for Matplotlib
    EXPECT_EQ(sizeOf(tiff), (std::vector<double>{2.0, 2.0}));
    auto u16 = [&](std::size_t at) { return unsigned(static_cast<unsigned char>(tiff[at])) | unsigned(static_cast<unsigned char>(tiff[at + 1])) << 8; };
    auto u32 = [&](std::size_t at) { return u16(at) | u16(at + 2) << 16; };
    std::map<unsigned, unsigned> tags;
    unsigned previous = 0;
    for (unsigned k = 0; k < u16(8); ++k) {
        const std::size_t e = 10 + 12 * k;
        EXPECT_GT(u16(e), previous);  // ascending, as TIFF requires
        previous = u16(e);
        tags[u16(e)] = u16(e + 2) == 3 && u32(e + 4) == 1 ? u16(e + 8) : u32(e + 8);
    }
    EXPECT_EQ(tags[259], 1u);  // uncompressed
    EXPECT_EQ(tags[262], 2u);  // RGB
    EXPECT_EQ(tags[277], 4u);  // with alpha
    EXPECT_EQ(tags[338], 2u);  // unassociated
    EXPECT_EQ(tags[279], 16u);
    EXPECT_EQ(tiff.substr(tags[273]), std::string(px.begin(), px.end()));
    EXPECT_EQ(u32(tags[282]), 3000000u);  // 300 dpi as 3000000 / 10000
    EXPECT_EQ(u32(tags[282] + 4), 10000u);
}

}  // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const std::unique_ptr<QGuiApplication> app = startBatchApplication(argv[0]);
    return RUN_ALL_TESTS();
}
