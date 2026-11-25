// The figure framework's QPainter backends (src/render/figure_qt.cpp): raster images of exactly the
// figure's size with the figure's colours, SVG documents with the texts kept as text and the
// figure's physical size, PDF pages of the figure's size, the frame loop's static and animated
// layers, files by extension. Runs on the offscreen platform (no display).

#include <array>
#include <cmath>
#include <tuple>
#include <filesystem>
#include <memory>
#include <regex>
#include <string>

#include <gtest/gtest.h>

#include <QGuiApplication>
#include <QImage>

#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "mdbridge/render/theme.hpp"

using namespace mdbridge::render;

namespace {

Figure demoFigure(int w, int h, double dpi) {
    Figure fig(w, h, dpi);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const Theme& th = lightTheme();
    fig.facecolor = toRgba(th.page);
    fig.grid.rows = 1;
    fig.grid.cols = 2;
    const std::size_t a = fig.addSubplot({0, 1, 0, 1});
    styleAxes(fig.axes[a], axesStyle(th, true, 11.0));
    Line l;
    l.x = {0.0, 1.0, 2.0, 3.0};
    l.y = {1.0, 3.0, 2.0, 5.0};
    l.color = toRgba(th.seriesColor(0));
    l.width = 1.6;
    l.zorder = 3;
    l.label = "RMSD";
    fig.axes[a].add(l);
    Line cursor;
    cursor.coords = Coords::XAxis;
    cursor.x = {1.5, 1.5};
    cursor.y = {0.0, 1.0};
    cursor.color = toRgba(th.ink);
    cursor.animated = true;
    fig.axes[a].add(cursor);
    fig.axes[a].setXLim(0.0, 3.0);
    fig.axes[a].xaxis.label.text = "Time (ns)";
    fig.axes[a].setTitle("RMSD (Å)", HAlign::Left).font.size = 11.0;
    Legend leg;
    leg.loc = LegendLoc::LowerRight;
    LegendEntry e;
    e.text = "RMSD";
    e.handle = l;
    leg.entries.push_back(e);
    styleLegend(leg, legendStyle(th, 11.0));
    fig.axes[a].legend = leg;
    // a heatmap with a colour bar
    const std::size_t b = fig.addSubplot({0, 1, 1, 2});
    styleAxes(fig.axes[b], axesStyle(th, false, 11.0));
    Image im;
    im.rows = 3;
    im.cols = 4;
    im.values = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, std::nan("")};
    im.cmap = th.sequential().withExtremes(toRgba(th.grid));
    im.vmin = 0;
    im.vmax = 10;
    im.extent = {-0.5, 3.5, 2.5, -0.5};
    fig.axes[b].add(im);
    fig.colorbar(b, 0, *im.cmap, 0.0, 10.0, 0.04, 0.02);
    fig.setSuptitle("Figure title", toRgba(th.ink), 11.0 * 1.3);
    Text& wm = fig.addText(0.997, 0.003, "Example data – layout preview");
    wm.ha = HAlign::Right;
    wm.va = VAlign::Bottom;
    wm.color = toRgba(th.muted);
    wm.font.size = 11.0 * 0.75;
    fig.draw(bundledFontSet());
    return fig;
}

std::array<int, 4> pixel(const QImage& img, int x, int y) {
    const QColor c = img.pixelColor(x, y);
    return {c.red(), c.green(), c.blue(), c.alpha()};
}

}  // namespace

TEST(FigureRaster, ExactSizeAndColours) {
    for (auto [w, h, dpi] : {std::tuple{960, 540, 60.0}, std::tuple{1000, 562, 77.0}, std::tuple{320, 180, 120.0}}) {
        const Figure fig = demoFigure(w, h, dpi);
        const QImage img = renderImage(buildDisplayList(fig, bundledFontSet()));
        EXPECT_EQ(img.width(), w);
        EXPECT_EQ(img.height(), h);
        EXPECT_EQ(img.format(), QImage::Format_RGBA8888);
        // the page colour in the corner, the axes surface inside the axes (away from lines and text)
        EXPECT_EQ(pixel(img, 1, h / 2), (std::array<int, 4>{0xf9, 0xf9, 0xf7, 255}));
        int inked = 0;
        for (int y = 0; y < h; y += 2) {
            for (int x = 0; x < w; x += 2) {
                if (img.pixelColor(x, y).value() < 100) ++inked;
            }
        }
        EXPECT_GT(inked, 10) << "text and lines drawn";
    }
}

TEST(FigureRaster, ImagesNearestAndMasked) {
    const Figure fig = demoFigure(960, 540, 60.0);
    const QImage img = renderImage(buildDisplayList(fig, bundledFontSet()));
    const Box b = fig.axesBox(1);
    // cell (row 0, col 0): value 0 = the ramp's first colour (the surface); last cell masked: grid colour
    const double cw = b.width() / 4, ch = b.height() / 3;
    const int x0 = static_cast<int>(b.x0 + cw / 2), y0 = static_cast<int>(540 - b.y1 + ch / 2);
    const Rgba c0 = lightTheme().sequential()(0.0);
    const auto p0 = pixel(img, x0, y0);
    EXPECT_NEAR(p0[0], c0.r * 255, 1.5);
    EXPECT_NEAR(p0[1], c0.g * 255, 1.5);
    EXPECT_NEAR(p0[2], c0.b * 255, 1.5);
    const int x3 = static_cast<int>(b.x0 + 3.5 * cw), y2 = static_cast<int>(540 - b.y1 + 2.5 * ch);
    EXPECT_EQ(pixel(img, x3, y2), (std::array<int, 4>{0xe1, 0xe0, 0xd9, 255}));
}

TEST(FigureRaster, LayersComposeToTheFrame) {
    // background (static) + the animated cursor = a frame; the background lacks the cursor
    const Figure fig = demoFigure(640, 360, 60.0);
    const QImage bg = renderImage(buildDisplayList(fig, bundledFontSet(), DisplayLayer::Static));
    const QImage frame = composeLayers(bg, buildDisplayList(fig, bundledFontSet(), DisplayLayer::Animated));
    ASSERT_EQ(frame.size(), bg.size());
    const Box b = fig.axesBox(0);
    const int x = static_cast<int>(std::lround(b.x0 + 0.5 * b.width()));
    const int y = static_cast<int>(360 - (b.y0 + 0.5 * b.height()));
    EXPECT_NE(pixel(frame, x, y), pixel(bg, x, y));
    int changed = 0;
    for (int yy = 0; yy < 360; ++yy) {
        for (int xx = 0; xx < 640; ++xx) {
            if (frame.pixel(xx, yy) != bg.pixel(xx, yy)) {
                ++changed;
                EXPECT_LE(std::abs(xx - x), 2) << "only the cursor changes";
            }
        }
    }
    EXPECT_GT(changed, 0);
}

TEST(FigureSvg, TextAsTextAndPhysicalSize) {
    const Figure fig = demoFigure(960, 540, 60.0);
    const QByteArray svg = renderSvg(buildDisplayList(fig, bundledFontSet()), QStringLiteral("preview"));
    const std::string s = svg.toStdString();
    EXPECT_NE(s.find("<svg"), std::string::npos);
    // 960 x 540 px at 60 dpi = 16 x 9 in = 406.4 x 228.6 mm
    std::smatch m;
    ASSERT_TRUE(std::regex_search(s, m, std::regex("width=\"([0-9.]+)mm\" height=\"([0-9.]+)mm\"")));
    EXPECT_NEAR(std::stod(m[1]), 406.4, 0.2);
    EXPECT_NEAR(std::stod(m[2]), 228.6, 0.2);
    EXPECT_NE(s.find("viewBox=\"0 0 960 540\""), std::string::npos);
    EXPECT_NE(s.find("<text"), std::string::npos);
    EXPECT_NE(s.find("DejaVu Sans"), std::string::npos);
    EXPECT_NE(s.find("Time (ns)"), std::string::npos);
    EXPECT_NE(s.find("Figure title"), std::string::npos);
    EXPECT_NE(s.find("RMSD (Å)"), std::string::npos);
}

TEST(FigurePdf, PageOfTheFigureSize) {
    const Figure fig = demoFigure(960, 540, 60.0);
    const QByteArray pdf = renderPdf(buildDisplayList(fig, bundledFontSet()));
    const std::string s = pdf.toStdString();
    ASSERT_EQ(s.rfind("%PDF", 0), 0u);
    // 16 x 9 in = 1152 x 648 pt
    std::smatch m;
    ASSERT_TRUE(std::regex_search(s, m, std::regex("/MediaBox \\[0 0 ([0-9.]+) ([0-9.]+)\\]")));
    EXPECT_NEAR(std::stod(m[1]), 1152.0, 0.5);
    EXPECT_NEAR(std::stod(m[2]), 648.0, 0.5);
    EXPECT_NE(s.find("/FontFile2"), std::string::npos) << "the bundled font embedded: text kept as text";
}

TEST(FigureFiles, ByExtension) {
    const Figure fig = demoFigure(320, 180, 60.0);
    const DisplayList dl = buildDisplayList(fig, bundledFontSet());
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mdbridge_figure_qt_test";
    std::filesystem::create_directories(dir);
    for (const char* ext : {"png", "svg", "pdf"}) {
        const std::filesystem::path p = dir / (std::string("figure.") + ext);
        writeFigure(dl, p);
        EXPECT_GT(std::filesystem::file_size(p), 100u) << ext;
    }
    const QImage png(QString::fromStdString((dir / "figure.png").string()));
    EXPECT_EQ(png.width(), 320);
    EXPECT_EQ(png.height(), 180);
    EXPECT_THROW(writeFigure(dl, dir / "figure.notaformat"), std::runtime_error);
    EXPECT_TRUE(canWriteFigureFormat("png"));
    EXPECT_TRUE(canWriteFigureFormat("eps"));  // render core's PostScript writer
    EXPECT_TRUE(canWriteFigureFormat("tiff"));
    EXPECT_FALSE(canWriteFigureFormat("notaformat"));
    std::filesystem::remove_all(dir);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const std::unique_ptr<QGuiApplication> app = startBatchApplication(argv[0]);
    return RUN_ALL_TESTS();
}
