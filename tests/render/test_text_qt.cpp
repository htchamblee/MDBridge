// The figure fonts in Qt (src/render/text_qt.cpp): the bundled faces from the program's resources,
// their registration with QFontDatabase, text drawn from the outlines, and the batch application on
// the "offscreen" platform (Linux: no display needed; CI runs this without X).

#include <memory>

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
    const QStringList families = render::registerBundledFonts();
    EXPECT_TRUE(families.contains("DejaVu Sans")) << families.join(", ").toStdString();
    EXPECT_TRUE(families.contains("DejaVu Sans Mono")) << families.join(", ").toStdString();
    EXPECT_EQ(render::registerBundledFonts(), families);  // once

    QImage image(10, 10, QImage::Format_RGBA8888_Premultiplied);
    image.setDotsPerMeterX(3937);  // 100 dpi
    image.setDotsPerMeterY(3937);
    const QFont font = render::figureFont(FontFace::MonoBold, 100.0 / 72.0 * 10.0, image);
    EXPECT_DOUBLE_EQ(font.pointSizeF(), 10.0 * 100.0 / static_cast<double>(image.logicalDpiY()));
    EXPECT_EQ(QFontInfo(font).family(), "DejaVu Sans Mono");
    EXPECT_TRUE(QFontInfo(font).bold());
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
    const render::TextLayout t = set.layout(FontFace::Sans, "RMSD (Å) Mg²⁺", 12, 150);
    QImage image(400, 80, QImage::Format_RGBA8888_Premultiplied);
    image.fill(Qt::transparent);
    QPainter painter(&image);
    render::drawText(painter, QPointF(10, 50), FontFace::Sans, t, Qt::black);
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
