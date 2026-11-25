#include "mdbridge/render/text_qt.hpp"

#include <mutex>
#include <string>

#include <QByteArray>
#include <QFile>
#include <QFontDatabase>
#include <QGuiApplication>
#include <QPaintDevice>
#include <QPainter>

#include "mdbridge/core/hash.hpp"

// the resources of qt_add_resources(... mdbridge_fonts ...) in a static library: initialised here
// so they are linked in whichever program uses the fonts (a second registration is ignored)
static void initFontResources() { Q_INIT_RESOURCE(mdbridge_fonts); }

namespace mdbridge::render {

namespace {

std::string resourceBytes(std::string_view file) {
    QFile f(QString::fromUtf8(kFontResourcePrefix.data(), static_cast<qsizetype>(kFontResourcePrefix.size())) +
            QString::fromUtf8(file.data(), static_cast<qsizetype>(file.size())));
    if (!f.open(QIODevice::ReadOnly)) {
        throw FontError("bundled font " + std::string(file) + " is missing from the program's resources");
    }
    const QByteArray bytes = f.readAll();
    return std::string(bytes.constData(), static_cast<std::size_t>(bytes.size()));
}

}  // namespace

const FontSet& bundledFontSet() {
    static const FontSet set = [] {
        initFontResources();
        return FontSet::load(resourceBytes);
    }();
    return set;
}

QStringList registerBundledFonts() {
    static std::mutex mutex;
    static QStringList families;
    const std::lock_guard lock(mutex);
    if (!families.isEmpty()) {
        return families;
    }
    const FontSet& set = bundledFontSet();
    QStringList found;
    for (const BundledFont& f : bundledFonts()) {
        const std::string& data = set.font(f.face).data();
        const int id = QFontDatabase::addApplicationFontFromData(
            QByteArray(data.data(), static_cast<qsizetype>(data.size())));
        const QStringList names = id < 0 ? QStringList() : QFontDatabase::applicationFontFamilies(id);
        if (names.isEmpty()) {
            throw FontError("Qt could not load the bundled font " + std::string(f.file));
        }
        for (const QString& name : names) {
            if (!found.contains(name)) {
                found.push_back(name);
            }
        }
    }
    families = found;
    return families;
}

QFont figureFont(FontFace face, double pixels, const QPaintDevice& device) {
    const BundledFont& f = bundledFont(face);
    QFont font(QString::fromUtf8(f.family.data(), static_cast<qsizetype>(f.family.size())));
    font.setWeight(f.weight >= 700 ? QFont::Bold : QFont::Normal);
    font.setStyle(f.oblique ? QFont::StyleOblique : QFont::StyleNormal);
    font.setFixedPitch(f.monospace);
    font.setPointSizeF(qtPointSize(pixels, device.logicalDpiY() > 0 ? device.logicalDpiY() : 96));
    font.setHintingPreference(QFont::PreferNoHinting);
    font.setStyleStrategy(QFont::StyleStrategy(QFont::NoFontMerging | QFont::PreferAntialias));
    font.setKerning(true);
    return font;
}

QPainterPath textPainterPath(FontFace face, const TextLayout& layout) {
    QPainterPath path;
    path.setFillRule(Qt::WindingFill);
    for (const PathCommand& c : textPath(bundledFontSet().font(face), layout)) {
        switch (c.kind) {
        case PathCommand::Kind::Move:
            path.moveTo(c.x, c.y);
            break;
        case PathCommand::Kind::Line:
            path.lineTo(c.x, c.y);
            break;
        case PathCommand::Kind::Quad:
            path.quadTo(c.cx, c.cy, c.x, c.y);
            break;
        case PathCommand::Kind::Close:
            path.closeSubpath();
            break;
        }
    }
    return path;
}

void drawText(QPainter& painter, QPointF origin, FontFace face, const TextLayout& layout, const QColor& color) {
    painter.save();
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.translate(origin);
    painter.fillPath(textPainterPath(face, layout), color);
    painter.restore();
}

QImage textSampleImage() {
    QImage image(720, 400, QImage::Format_RGBA8888_Premultiplied);
    image.fill(Qt::white);
    QPainter painter(&image);
    const FontSet& set = bundledFontSet();
    const char* const sample = "MDBridge: RMSD (\u00c5) Mg\u00b2\u207a \u22120.5 \u00d710\u00b3 \u0394t fi AV";
    double y = 28;
    for (const BundledFont& f : bundledFonts()) {
        for (const double points : {7.0, 10.0, 13.0}) {
            const TextLayout t = set.layout(f.face, sample, points, 100);
            drawText(painter, QPointF(40.25, y), f.face, t, QColor(0x1f, 0x2a, 0x44));
            y += t.pixelSize * 1.3;
        }
    }
    // a y-axis label (rotated) and coloured, translucent text
    painter.save();
    painter.translate(18, 300);
    painter.rotate(-90);
    drawText(painter, QPointF(0, 0), FontFace::Sans, set.layout(FontFace::Sans, "Time (ns)", 11, 100), QColor(0xd6, 0x27, 0x28));
    painter.restore();
    drawText(painter, QPointF(500.5, 380.7), FontFace::SansBold, set.layout(FontFace::SansBold, "frame 120", 12, 150),
             QColor(31, 119, 180, 160));
    painter.end();
    return image;
}

std::string imagePixelsSha256(const QImage& image) {
    const QImage rgba = image.convertToFormat(QImage::Format_RGBA8888_Premultiplied);
    std::string bytes;
    bytes.reserve(static_cast<std::size_t>(rgba.width()) * 4 * static_cast<std::size_t>(rgba.height()));
    for (int row = 0; row < rgba.height(); ++row) {
        bytes.append(reinterpret_cast<const char*>(rgba.constScanLine(row)), static_cast<std::size_t>(rgba.width()) * 4);
    }
    return core::sha256Hex(bytes);
}

std::unique_ptr<QGuiApplication> startBatchApplication(const char* program) {
    if (const auto platform = batchQtPlatform()) {
        qputenv("QT_QPA_PLATFORM", QByteArray(platform->data(), static_cast<qsizetype>(platform->size())));
    }
    // Qt keeps references to argc / argv for the application's lifetime
    static int argc = 1;
    static std::string name;
    static char* argv[] = {nullptr, nullptr};
    name = program ? program : "mdbridge";
    argv[0] = name.data();
    return std::make_unique<QGuiApplication>(argc, argv);
}

}  // namespace mdbridge::render
