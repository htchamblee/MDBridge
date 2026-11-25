#include "software_view.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>

#include <QPainter>
#include <QPainterPath>
#include <QRadialGradient>
#include <QThreadPool>

#include "scene.hpp"

namespace mdbridge::render {

namespace {

constexpr int kFrameIntervalMs = 33;  // ~30 frames per second is plenty for a software view

QColor rgbColour(const std::array<float, 3>& rgb) {
    return QColor::fromRgbF(std::clamp(rgb[0], 0.0f, 1.0f), std::clamp(rgb[1], 0.0f, 1.0f), std::clamp(rgb[2], 0.0f, 1.0f));
}

// text that reads on `background`
QColor inkOn(const QColor& background) {
    return background.lightnessF() > 0.5 ? QColor(70, 76, 86) : QColor(170, 178, 189);
}

// The molecule over the background, `width` x `height` device pixels.
QImage renderMolecule(const MoleculeFrame& frame, const ViewCamera& camera, int width, int height) {
    const molecule::SoftwareRenderer renderer;
    const molecule::RgbaImage img =
        renderer.render(frame.items, cameraViewport(frame, camera, width, height), frame.shading, camera.quality);
    QImage layer(static_cast<int>(img.width), static_cast<int>(img.height), QImage::Format_RGBA8888);
    for (int y = 0; y < layer.height(); ++y) {
        std::copy_n(img.rgba.data() + static_cast<std::size_t>(y) * img.width * 4, img.width * 4, layer.scanLine(y));
    }
    QImage out(width, height, QImage::Format_RGB32);
    out.fill(rgbColour(frame.background));
    QPainter painter(&out);
    painter.drawImage(QRectF(0, 0, width, height), layer);
    return out;
}

}  // namespace

struct SoftwareWindow::Receiver {
    std::mutex mutex;
    SoftwareWindow* window = nullptr;
};

SoftwareWindow::SoftwareWindow(QString note) : m_note(std::move(note)), m_receiver(std::make_shared<Receiver>()) {
    m_receiver->window = this;
    m_clock.start();
    m_timer.start(kFrameIntervalMs, this);
}

SoftwareWindow::~SoftwareWindow() {
    const std::lock_guard<std::mutex> lock(m_receiver->mutex);
    m_receiver->window = nullptr;  // a render still running drops its image
}

QImage SoftwareWindow::render(QSize size, qreal devicePixelRatio) const {
    QImage image(size * devicePixelRatio, QImage::Format_ARGB32_Premultiplied);
    image.setDevicePixelRatio(devicePixelRatio);
    QPainter painter(&image);
    draw(painter, size);
    return image;
}

void SoftwareWindow::setMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera) {
    m_frame = std::move(frame);
    m_camera = camera;
    ++m_generation;
    m_problem.clear();
    if (m_frame) {
        m_timer.stop();
        scheduleRender();
    } else {
        m_image = QImage();
        if (!m_timer.isActive()) {
            m_timer.start(kFrameIntervalMs, this);
        }
    }
    update();
}

bool SoftwareWindow::moleculeShown() const {
    return m_frame && !m_image.isNull() && m_problem.isEmpty() && m_imageGeneration == m_generation;
}

void SoftwareWindow::scheduleRender() {
    if (!m_frame) {
        return;
    }
    if (m_busy) {
        m_pending = true;  // the newest request waits; older ones are forgotten
        return;
    }
    startRender();
}

void SoftwareWindow::startRender() {
    m_pending = false;
    if (!m_frame) {
        return;
    }
    const qreal dpr = devicePixelRatio();
    const int width = static_cast<int>(std::lround(this->width() * dpr));
    const int height = static_cast<int>(std::lround(this->height() * dpr));
    if (width < 2 || height < 2) {
        return;  // not laid out yet: the resize renders
    }
    m_busy = true;
    auto job = [frame = m_frame, camera = m_camera, width, height, dpr, generation = m_generation,
                receiver = m_receiver]() {
        QImage image;
        QString problem;
        try {
            image = renderMolecule(*frame, camera, width, height);
            image.setDevicePixelRatio(dpr);
        } catch (const std::exception& e) {
            problem = QString::fromUtf8(e.what());
        }
        const std::lock_guard<std::mutex> lock(receiver->mutex);
        if (SoftwareWindow* w = receiver->window) {
            QMetaObject::invokeMethod(
                w, [w, generation, image, problem] { w->rendered(generation, image, problem); }, Qt::QueuedConnection);
        }
    };
    QThreadPool::globalInstance()->start(std::function<void()>(std::move(job)));
}

void SoftwareWindow::rendered(std::uint64_t generation, const QImage& image, const QString& problem) {
    m_busy = false;
    if (m_frame && generation >= m_imageGeneration) {
        m_imageGeneration = generation;
        m_problem = problem;
        if (!image.isNull()) {
            m_image = image;
        }
    }
    if (m_pending || (m_frame && generation != m_generation)) {
        startRender();
    }
    update();
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

void SoftwareWindow::resizeEvent(QResizeEvent* event) {
    QRasterWindow::resizeEvent(event);
    if (m_frame) {
        ++m_generation;
        scheduleRender();
    }
}

void SoftwareWindow::draw(QPainter& painter, QSizeF size) const {
    if (!m_frame) {
        drawPlaceholder(painter, size);
        return;
    }
    const QColor background = rgbColour(m_frame->background);
    painter.fillRect(QRectF(QPointF(0, 0), size), background);
    const QColor ink = inkOn(background);
    painter.setPen(ink);
    const QRectF text(QPointF(10, 0), size - QSizeF(20, 10));
    if (!m_image.isNull()) {
        painter.setRenderHint(QPainter::SmoothPixmapTransform);
        painter.drawImage(QRectF(QPointF(0, 0), size), m_image);  // scaled until the next image is ready
    } else if (m_problem.isEmpty()) {
        painter.drawText(QRectF(QPointF(0, 0), size), Qt::AlignCenter, tr("Drawing the structure..."));
    }
    if (!m_problem.isEmpty()) {
        painter.drawText(QRectF(QPointF(0, 0), size), Qt::AlignCenter | Qt::TextWordWrap,
                         tr("The structure could not be drawn: %1").arg(m_problem));
    }
    painter.drawText(text, Qt::AlignLeft | Qt::AlignBottom, tr("Software view: %1").arg(m_note));
}

void SoftwareWindow::drawPlaceholder(QPainter& painter, QSizeF size) const {
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

bool SoftwareView::showMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera) {
    auto* w = static_cast<SoftwareWindow*>(window());
    if (!w) {
        return false;
    }
    w->setMolecule(std::move(frame), camera);
    return true;
}

bool SoftwareView::moleculeShown() const {
    const auto* w = static_cast<const SoftwareWindow*>(window());
    return w && w->moleculeShown();
}

QImage SoftwareView::grabFrame() {
    auto* w = static_cast<SoftwareWindow*>(window());
    if (!w || w->size().isEmpty()) {
        return {};
    }
    return w->render(w->size(), w->devicePixelRatio());
}

}  // namespace mdbridge::render
