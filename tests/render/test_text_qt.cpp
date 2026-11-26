// The bundled fonts in Qt (src/render/text_qt.cpp): the bundled faces from the program's resources,
// their registration with QFontDatabase, the interface's fonts, the font runs of vector text, text
// drawn from the outlines, and the batch application on the "offscreen" platform (Linux: no display
// needed; CI runs this without X).

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <QFontInfo>
#include <QGuiApplication>
#include <QImage>
#include <QPainter>

#include "mdbridge/render/text_qt.hpp"

using namespace mdbridge;
using render::FontFace;

TEST(TextQt, BundledFontsFromResources) {
    const render::FontSet& set = render::bundledFontSet();
    for (const render::BundledFont& f : render::bundledFonts()) {
        EXPECT_EQ(set.font(f.face).family(), f.family);
        EXPECT_EQ(set.font(f.face).subfamily(), f.style);
    }
    EXPECT_EQ(set.fallback().family(), render::fallbackFont().family);
    const QStringList families = render::registerBundledFonts();
    for (const char* family : {"IBM Plex Serif", "IBM Plex Math", "IBM Plex Mono"}) {
        EXPECT_TRUE(families.contains(family)) << families.join(", ").toStdString();
    }
    EXPECT_EQ(render::registerBundledFonts(), families);  // once

    QImage image(10, 10, QImage::Format_RGBA8888_Premultiplied);
    image.setDotsPerMeterX(3937);  // 100 dpi
    image.setDotsPerMeterY(3937);
    const QFont font = render::figureFont(FontFace::BoldItalic, 100.0 / 72.0 * 10.0, image);
    EXPECT_DOUBLE_EQ(font.pointSizeF(), 10.0 * 100.0 / static_cast<double>(image.logicalDpiY()));
    EXPECT_EQ(QFontInfo(font).family(), "IBM Plex Serif");
    EXPECT_TRUE(QFontInfo(font).bold());
    EXPECT_TRUE(QFontInfo(font).italic());
    const QFont fallback = render::figureFont(FontFace::BoldItalic, 100.0 / 72.0 * 10.0, image, true);
    EXPECT_EQ(QFontInfo(fallback).family(), "IBM Plex Math");
    EXPECT_FALSE(QFontInfo(fallback).bold());  // its one face, not a synthesized bold

    // the interface: Plex Serif, Plex Mono for fixed-width text, Plex Math next in both
    const QFont ui = render::interfaceFont(9);
    EXPECT_EQ(QFontInfo(ui).family(), "IBM Plex Serif");
    EXPECT_EQ(ui.families(), QStringList({"IBM Plex Serif", "IBM Plex Math"}));
    EXPECT_DOUBLE_EQ(ui.pointSizeF(), 9.0);
    const QFont fixed = render::interfaceFixedFont(9);
    EXPECT_EQ(QFontInfo(fixed).family(), "IBM Plex Mono");
    EXPECT_TRUE(QFontInfo(fixed).fixedPitch());
    EXPECT_EQ(fixed.families(), QStringList({"IBM Plex Mono", "IBM Plex Math"}));
}

TEST(TextQt, FontRunsOfVectorText) {
    const render::FontSet& set = render::bundledFontSet();
    const std::string text = "Mg\u00b2\u207a ions \u207b";
    const render::TextLayout t = set.layout(FontFace::Regular, text, 10, 100);
    const std::vector<render::TextRun> runs = render::fontRuns(text, t);
    ASSERT_EQ(runs.size(), 4u);  // "Mg²", "⁺", " ions ", "⁻"
    EXPECT_FALSE(runs[0].fallback);
    EXPECT_TRUE(runs[1].fallback);
    EXPECT_FALSE(runs[2].fallback);
    EXPECT_TRUE(runs[3].fallback);
    EXPECT_EQ(render::runText(text, runs[0]), QStringLiteral("Mg\u00b2"));
    EXPECT_EQ(render::runText(text, runs[1]), QStringLiteral("\u207a"));
    EXPECT_EQ(render::runText(text, runs[2]), QStringLiteral(" ions "));
    EXPECT_EQ(render::runText(text, runs[3]), QStringLiteral("\u207b"));
    EXPECT_DOUBLE_EQ(runs[0].x, 0.0);
    EXPECT_DOUBLE_EQ(runs[1].x, static_cast<double>(t.glyphs[3].x) / 64.0);
    EXPECT_TRUE(render::fontRuns("", set.layout(FontFace::Regular, "", 10, 100)).empty());
}

TEST(TextQt, BatchPlatform) {
#if defined(_WIN32) || defined(__APPLE__)
    EXPECT_NE(QGuiApplication::platformName(), "");
#else
    EXPECT_EQ(QGuiApplication::platformName(), "offscreen");
#endif
}

TEST(TextQt, DrawsTheLayoutsOutlines) {
    const render::FontSet& set = render::bundledFontSet();
    const render::TextLayout t = set.layout(FontFace::Regular, "RMSD (Å) Mg²⁺", 12, 150);
    QImage image(400, 80, QImage::Format_RGBA8888_Premultiplied);
    image.fill(Qt::transparent);
    QPainter painter(&image);
    render::drawText(painter, QPointF(10, 50), FontFace::Regular, t, Qt::black);
    painter.end();
    // ink only inside the layout's ink box (plus one pixel of antialiasing), and plenty of it
    int inked = 0;
    for (int y = 0; y < image.height(); ++y) {
        for (int x = 0; x < image.width(); ++x) {
            if (qAlpha(image.pixel(x, y)) == 0) {
                continue;
            }
            ++inked;
            EXPECT_GE(x, 10 + t.inkXMin / 64.0 - 1);
            EXPECT_LE(x, 10 + t.inkXMax / 64.0 + 1);
            EXPECT_GE(y, 50 - t.inkYMax / 64.0 - 1);
            EXPECT_LE(y, 50 - t.inkYMin / 64.0 + 1);
        }
    }
    EXPECT_GT(inked, 600);
}

TEST(TextQt, SampleImageIsReproducible) {
    const QImage a = render::textSampleImage();
    const QImage b = render::textSampleImage();
    EXPECT_EQ(a.size(), QSize(720, 400));
    EXPECT_EQ(render::imagePixelsSha256(a), render::imagePixelsSha256(b));
    EXPECT_EQ(render::imagePixelsSha256(a).size(), 64u);
    int dark = 0;
    for (int y = 0; y < a.height(); ++y) {
        for (int x = 0; x < a.width(); ++x) {
            dark += qGray(a.pixel(x, y)) < 128;
        }
    }
    EXPECT_GT(dark, 5000);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const std::unique_ptr<QGuiApplication> app = render::startBatchApplication(argv[0]);
    return RUN_ALL_TESTS();
}
