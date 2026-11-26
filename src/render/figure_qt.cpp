#include "mdbridge/render/figure_qt.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <QBuffer>
#include <QColor>
#include <QFile>
#include <QImageWriter>
#include <QPageSize>
#include <QPainter>
#include <QPainterPath>
#include <QPdfWriter>
#include <QPen>
#include <QSvgGenerator>
#include <QTransform>

#include "mdbridge/render/figure_formats.hpp"
#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::render {

namespace {

QColor qcolor(const Rgba& c) {
    QColor q;
    q.setRgbF(static_cast<float>(std::clamp(c.r, 0.0, 1.0)), static_cast<float>(std::clamp(c.g, 0.0, 1.0)),
              static_cast<float>(std::clamp(c.b, 0.0, 1.0)), static_cast<float>(std::clamp(c.a, 0.0, 1.0)));
    return q;
}

QPen qpen(const DisplayStroke& s) {
    QPen pen(qcolor(s.color));
    pen.setWidthF(s.width);
    pen.setCapStyle(s.cap == CapStyle::Butt ? Qt::FlatCap : s.cap == CapStyle::Round ? Qt::RoundCap : Qt::SquareCap);
    pen.setJoinStyle(s.roundJoin ? Qt::RoundJoin : Qt::MiterJoin);
    if (!s.dashes.empty() && s.width > 0) {
        // QPen's dash pattern is in units of the pen width; an odd pattern repeats (as Agg does)
        QList<qreal> pattern;
        for (int rep = 0; rep < (s.dashes.size() % 2 ? 2 : 1); ++rep) {
            for (double d : s.dashes) pattern.push_back(std::max(d, 1e-6) / s.width);
        }
        pen.setDashPattern(pattern);
        pen.setDashOffset(s.dashOffset / s.width);
    }
    return pen;
}

QRectF qrect(const Box& b) { return QRectF(QPointF(b.x0, b.y0), QPointF(b.x1, b.y1)); }

struct Painter {
    QPainter& p;
    TextOutput text;

    void clip(const std::optional<Box>& c) {
        if (c) {
            p.setClipRect(qrect(*c));
        } else {
            p.setClipping(false);
        }
    }

    void operator()(const FillRectItem& r) {
        p.setRenderHint(QPainter::Antialiasing, false);
        p.fillRect(qrect(r.rect), qcolor(r.color));
        p.setRenderHint(QPainter::Antialiasing, true);
    }

    void operator()(const PathItem& item) {
        QPainterPath path;
        path.setFillRule(Qt::WindingFill);
        for (const auto& sub : item.subpaths) {
            if (sub.empty()) continue;
            path.moveTo(sub[0].x, sub[0].y);
            for (std::size_t i = 1; i < sub.size(); ++i) path.lineTo(sub[i].x, sub[i].y);
            if (item.closed) path.closeSubpath();
        }
        if (item.fill) p.fillPath(path, qcolor(*item.fill));
        if (item.stroke && item.stroke->width > 0) p.strokePath(path, qpen(*item.stroke));
    }

    void operator()(const MarkerItem& m) {
        for (std::size_t i = 0; i < m.centers.size(); ++i) {
            const double size = m.sizes.empty() ? 0.0 : m.sizes[i % m.sizes.size()];
            const QPointF c(m.centers[i].x, m.centers[i].y);
            if (m.edges.empty() || m.edgeWidth <= 0) {
                p.setPen(Qt::NoPen);
            } else {
                QPen pen(qcolor(m.edges[i % m.edges.size()]));
                pen.setWidthF(m.edgeWidth);
                pen.setJoinStyle(Qt::MiterJoin);
                p.setPen(pen);
            }
            if (m.faces.empty()) {
                p.setBrush(Qt::NoBrush);
            } else {
                p.setBrush(qcolor(m.faces[i % m.faces.size()]));
            }
            if (m.shape == Marker::Square) {
                p.drawRect(QRectF(c.x() - size / 2, c.y() - size / 2, size, size));
            } else {
                p.drawEllipse(c, size / 2, size / 2);
            }
        }
    }

    void operator()(const RoundRectItem& r) {
        p.setPen(r.stroke && r.stroke->width > 0 ? qpen(*r.stroke) : QPen(Qt::NoPen));
        p.setBrush(r.fill ? QBrush(qcolor(*r.fill)) : QBrush(Qt::NoBrush));
        p.drawRoundedRect(qrect(r.rect), r.radius, r.radius);
    }

    void operator()(const TextItem& t) {
        p.save();
        p.translate(t.origin.x, t.origin.y);
        if (t.rotation != 0.0) p.rotate(-t.rotation);
        if (text == TextOutput::Outlines) {
            drawText(p, QPointF(0, 0), t.face, t.layout, qcolor(t.color));
        } else {
            // one QFont per run: the figure face, and the fallback font for the glyphs it lacks
            p.setPen(qcolor(t.color));
            for (const TextRun& run : fontRuns(t.text, t.layout)) {
                p.setFont(figureFont(t.face, t.pixelSize, *p.device(), run.fallback));
                p.drawText(QPointF(run.x, 0), runText(t.text, run));
            }
        }
        p.restore();
    }

    void operator()(const ImageItem& im) {
        if (im.width <= 0 || im.height <= 0) return;
        const QImage img(im.rgba.data(), im.width, im.height, 4 * im.width, QImage::Format_RGBA8888);
        p.save();
        p.setRenderHint(QPainter::Antialiasing, false);
        p.setRenderHint(QPainter::SmoothPixmapTransform, im.smooth);
        const double sx = (im.corner1.x - im.corner0.x) / im.width, sy = (im.corner1.y - im.corner0.y) / im.height;
        p.setTransform(QTransform(sx, 0, 0, sy, im.corner0.x, im.corner0.y), true);
        p.drawImage(QRectF(0, 0, im.width, im.height), img);
        p.restore();
    }
};

}  // namespace

void paintDisplayList(QPainter& painter, const DisplayList& list, TextOutput text) {
    painter.save();
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);
    Painter p{painter, text};
    for (const DisplayItem& item : list.items) {
        p.clip(item.clip);
        std::visit(p, item.item);
    }
    painter.restore();
}

QImage renderImage(const DisplayList& list) {
    QImage img(list.width, list.height, QImage::Format_ARGB32_Premultiplied);
    img.fill(Qt::transparent);
    {
        QPainter painter(&img);
        paintDisplayList(painter, list, TextOutput::Outlines);
    }
    return img.convertToFormat(QImage::Format_RGBA8888);
}

QImage composeLayers(const QImage& background, const DisplayList& layer) {
    QImage img = background.convertToFormat(QImage::Format_ARGB32_Premultiplied);
    {
        QPainter painter(&img);
        paintDisplayList(painter, layer, TextOutput::Outlines);
    }
    return img.convertToFormat(QImage::Format_RGBA8888);
}

QByteArray renderSvg(const DisplayList& list, const QString& title) {
    registerBundledFonts();
    QByteArray bytes;
    QBuffer buffer(&bytes);
    buffer.open(QIODevice::WriteOnly);
    QSvgGenerator gen;
    gen.setOutputDevice(&buffer);
    gen.setSize(QSize(list.width, list.height));
    gen.setViewBox(QRect(0, 0, list.width, list.height));
    gen.setResolution(std::max(1, static_cast<int>(std::lround(list.dpi))));
    if (!title.isEmpty()) gen.setTitle(title);
    gen.setDescription(QStringLiteral("MDBridge"));
    {
        QPainter painter(&gen);
        paintDisplayList(painter, list, TextOutput::Fonts);
    }
    return bytes;
}

QByteArray renderPdf(const DisplayList& list, const QString& title) {
    registerBundledFonts();
    QByteArray bytes;
    QBuffer buffer(&bytes);
    buffer.open(QIODevice::WriteOnly);
    {
        QPdfWriter pdf(&buffer);
        pdf.setCreator(QStringLiteral("MDBridge"));
        if (!title.isEmpty()) pdf.setTitle(title);
        // the page is the figure: width / dpi inches, no margins; the painter works in figure pixels
        const QSizeF points(list.width / list.dpi * 72.0, list.height / list.dpi * 72.0);
        pdf.setPageSize(QPageSize(points, QPageSize::Point, QString(), QPageSize::ExactMatch));
        pdf.setPageMargins(QMarginsF(0, 0, 0, 0));
        constexpr int kResolution = 1200;
        pdf.setResolution(kResolution);
        QPainter painter(&pdf);
        painter.scale(kResolution / list.dpi, kResolution / list.dpi);
        paintDisplayList(painter, list, TextOutput::Fonts);
    }
    return bytes;
}

namespace {
std::string lowerExtension(const std::filesystem::path& path) {
    std::string ext = path.extension().string();
    if (!ext.empty() && ext[0] == '.') ext.erase(0, 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return ext;
}
}  // namespace

bool canWriteFigureFormat(const std::string& format) {
    if (format == "png" || format == "svg" || format == "pdf" || format == "eps" || format == "ps" || format == "tif" ||
        format == "tiff") {
        return true;
    }
    for (const QByteArray& f : QImageWriter::supportedImageFormats()) {
        if (f.toStdString() == format) return true;
    }
    return false;
}

namespace {
void writeBytes(const std::filesystem::path& path, const char* data, qsizetype size) {
    QFile out(QString::fromStdWString(path.wstring()));
    if (!out.open(QIODevice::WriteOnly) || out.write(data, size) != size) {
        throw std::runtime_error("cannot write " + path.string());
    }
}
}  // namespace

void writeFigure(const DisplayList& list, const std::filesystem::path& path, const std::string& format) {
    const std::string ext = format.empty() ? lowerExtension(path) : format;
    if (ext == "svg" || ext == "pdf") {
        const QByteArray bytes = ext == "svg" ? renderSvg(list) : renderPdf(list);
        writeBytes(path, bytes.constData(), bytes.size());
        return;
    }
    if (ext == "eps" || ext == "ps") {
        const std::string ps = renderPostScript(list, bundledFontSet(), ext == "eps");
        writeBytes(path, ps.data(), static_cast<qsizetype>(ps.size()));
        return;
    }
    if (!canWriteFigureFormat(ext)) {
        throw std::runtime_error("cannot write " + path.string() + ": unknown image format '" + ext + "'");
    }
    QImage img = renderImage(list);
    if (ext == "tif" || ext == "tiff") {
        const std::string tiff = encodeTiff(img.width(), img.height(),
                                            std::span<const std::uint8_t>(img.constBits(), static_cast<std::size_t>(img.sizeInBytes())),
                                            list.dpi);
        writeBytes(path, tiff.data(), static_cast<qsizetype>(tiff.size()));
        return;
    }
    if (ext == "jpg" || ext == "jpeg") {
        // no alpha in JPEG: the figure over white, as Matplotlib's savefig does
        QImage flat(img.size(), QImage::Format_RGB32);
        flat.fill(Qt::white);
        QPainter(&flat).drawImage(0, 0, img);
        img = flat;
    }
    const int dotsPerMeter = static_cast<int>(std::lround(list.dpi / 0.0254));
    img.setDotsPerMeterX(dotsPerMeter);
    img.setDotsPerMeterY(dotsPerMeter);
    QImageWriter writer(QString::fromStdWString(path.wstring()), QByteArray::fromStdString(ext));
    if (!writer.write(img)) {
        throw std::runtime_error("cannot write " + path.string() + ": " + writer.errorString().toStdString());
    }
}

void writeFigure(const Figure& figure, const std::filesystem::path& path, const std::string& format) {
    writeFigure(buildDisplayList(figure, bundledFontSet()), path, format);
}

}  // namespace mdbridge::render
