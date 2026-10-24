#include "software_view.hpp"

#include <algorithm>
#include <cmath>

#include <QPainter>
#include <QPainterPath>
#include <QRadialGradient>

#include "scene.hpp"

namespace mdbridge::render {

namespace {
constexpr int kFrameIntervalMs = 33;  // ~30 frames per second is plenty for a software view
}

SoftwareWindow::SoftwareWindow(QString note) : m_note(std::move(note)) {
    m_clock.start();
    m_timer.start(kFrameIntervalMs, this);
}

QImage SoftwareWindow::render(QSize size, qreal devicePixelRatio) const {
    QImage image(size * devicePixelRatio, QImage::Format_ARGB32_Premultiplied);
    image.setDevicePixelRatio(devicePixelRatio);
    QPainter painter(&image);
    draw(painter, size);
    return image;
}

void SoftwareWindow::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    draw(painter, size());
    emit frameRendered();
}

void SoftwareWindow::timerEvent(QTimerEvent* event) {
    if (event->timerId() == m_timer.timerId() && isExposed()) {
        update();
    }
}

void SoftwareWindow::draw(QPainter& painter, QSizeF size) const {
    using namespace scene;
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(QRectF(QPointF(0, 0), size),
                     QColor::fromRgbF(kBackground[0], kBackground[1], kBackground[2]));

    // the corners in pixels, rotated and scaled as the Vulkan view's vertex shader does
    const double angle = angleAt(m_clock.elapsed());
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double aspect = size.height() > 0 ? size.width() / size.height() : 1.0;
    QPointF corners[3];
    for (int i = 0; i < 3; ++i) {
        const QPointF p = kCorners[i];
        const double x = (c * p.x() - s * p.y()) / std::max(aspect, 1e-3);
        const double y = s * p.x() + c * p.y();
        corners[i] = QPointF((x + 1.0) / 2.0 * size.width(), (y + 1.0) / 2.0 * size.height());
    }
    QPainterPath triangle;
    triangle.addPolygon(QPolygonF({corners[0], corners[1], corners[2], corners[0]}));

    // per-corner colours fading towards the middle, close to the Vulkan view's interpolation
    painter.fillPath(triangle, QColor(130, 140, 150));
    const double edge = std::hypot(corners[1].x() - corners[0].x(), corners[1].y() - corners[0].y());
    painter.setClipPath(triangle);
    for (int i = 0; i < 3; ++i) {
        QRadialGradient glow(corners[i], edge);
        QColor transparent = kColours[i];
        transparent.setAlpha(0);
        glow.setColorAt(0.0, kColours[i]);
        glow.setColorAt(1.0, transparent);
        painter.fillPath(triangle, glow);
    }
    painter.setClipping(false);

    painter.setPen(QColor(170, 178, 189));
    painter.drawText(QRectF(QPointF(10, 0), size - QSizeF(20, 10)), Qt::AlignLeft | Qt::AlignBottom,
                     tr("Software view: %1").arg(m_note));
}

SoftwareView::SoftwareView(const QString& reason)
    : View(new SoftwareWindow(reason), Renderer::Software, tr("Software view (%1)").arg(reason)) {
    connect(static_cast<SoftwareWindow*>(window()), &SoftwareWindow::frameRendered, this, [this] { countFrame(); });
}

QImage SoftwareView::grabFrame() {
    auto* w = static_cast<SoftwareWindow*>(window());
    if (!w || w->size().isEmpty()) {
        return {};
    }
    return w->render(w->size(), w->devicePixelRatio());
}

}  // namespace mdbridge::render
