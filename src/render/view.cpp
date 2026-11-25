#include "mdbridge/render/view.hpp"

#include <algorithm>
#include <cmath>

#include <QMouseEvent>
#include <QTimer>
#include <QWheelEvent>

#include "software_view.hpp"
#include "vulkan_view.hpp"

namespace mdbridge::render {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kSettleMs = 250;  // fast quality until the wheel has been still this long

using Mat = std::array<molecule::Vec3, 3>;

Mat multiply(const Mat& a, const Mat& b) {
    Mat out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    return out;
}

// Gram-Schmidt on the rows: many small turns keep the matrix a rotation
Mat orthonormal(Mat m) {
    auto dot = [](const molecule::Vec3& a, const molecule::Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
    auto normalise = [&](molecule::Vec3& v) {
        const double n = std::sqrt(dot(v, v));
        if (n > 0) {
            for (double& c : v) c /= n;
        }
    };
    normalise(m[0]);
    const double d = dot(m[1], m[0]);
    for (std::size_t k = 0; k < 3; ++k) m[1][k] -= d * m[0][k];
    normalise(m[1]);
    m[2] = {m[0][1] * m[1][2] - m[0][2] * m[1][1], m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0]};
    return m;
}

}  // namespace

ViewCamera frameCamera(const MoleculeFrame& frame) {
    ViewCamera camera;
    camera.rotation = frame.rotation;
    camera.quality = frame.quality;
    return camera;
}

molecule::Viewport cameraViewport(const MoleculeFrame& frame, const ViewCamera& camera, int width, int height) {
    molecule::Viewport vp;
    vp.rotation = camera.rotation;
    vp.zrange = frame.zrange;
    vp.width = std::max(width, 1);
    vp.height = std::max(height, 1);
    const double zoom = camera.zoom > 0 && std::isfinite(camera.zoom) ? camera.zoom : 1.0;
    const double dx = (frame.xlim[1] - frame.xlim[0]) / zoom, dy = (frame.ylim[1] - frame.ylim[0]) / zoom;
    const double cx = 0.5 * (frame.xlim[0] + frame.xlim[1]), cy = 0.5 * (frame.ylim[0] + frame.ylim[1]);
    const double target = static_cast<double>(vp.width) / static_cast<double>(vp.height);
    double hx = 0.5 * dx, hy = 0.5 * dy;
    if (dy > 0 && dx / dy < target) {
        hx = hy * target;  // window wider than the data: widen x
    } else if (dx > 0) {
        hy = hx / target;
    }
    vp.xlim = {cx - hx, cx + hx};
    vp.ylim = {cy - hy, cy + hy};
    return vp;
}

molecule::Viewport fittedViewport(const MoleculeFrame& frame, int width, int height) {
    return cameraViewport(frame, frameCamera(frame), width, height);
}

ViewCamera rotatedCamera(const ViewCamera& camera, double dx, double dy, double height) {
    const double perPixel = kPi / std::max(height, 1.0);
    const double a = dx * perPixel;  // about the screen's vertical axis: the front (+z) follows to the right
    const double b = dy * perPixel;  // about the horizontal axis: the front follows down (view y is up)
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b);
    const Mat ry{{{ca, 0.0, sa}, {0.0, 1.0, 0.0}, {-sa, 0.0, ca}}};
    const Mat rx{{{1.0, 0.0, 0.0}, {0.0, cb, -sb}, {0.0, sb, cb}}};
    ViewCamera out = camera;
    out.rotation = orthonormal(multiply(rx, multiply(ry, camera.rotation)));
    return out;
}

ViewCamera zoomedCamera(const ViewCamera& camera, double steps) {
    ViewCamera out = camera;
    out.zoom = std::clamp(camera.zoom * std::pow(1.15, steps), 0.05, 50.0);
    return out;
}

std::optional<Renderer> rendererFromName(std::string_view name) {
    if (name == "auto") {
        return Renderer::Auto;
    }
    if (name == "vulkan") {
        return Renderer::Vulkan;
    }
    if (name == "software") {
        return Renderer::Software;
    }
    return std::nullopt;
}

const char* rendererName(Renderer renderer) {
    switch (renderer) {
    case Renderer::Auto:
        return "auto";
    case Renderer::Vulkan:
        return "vulkan";
    case Renderer::Software:
        return "software";
    }
    return "auto";
}

View::View(QWindow* window, Renderer renderer, QString description)
    : m_window(window), m_renderer(renderer), m_description(std::move(description)) {
    m_settle = new QTimer(this);
    m_settle->setSingleShot(true);
    m_settle->setInterval(kSettleMs);
    connect(m_settle, &QTimer::timeout, this, &View::settle);
    if (window) {
        window->installEventFilter(this);
    }
}

bool View::setMoleculeFrame(std::shared_ptr<const MoleculeFrame> frame) {
    m_dragging = false;
    m_settle->stop();
    m_camera = frame ? frameCamera(*frame) : ViewCamera();
    m_frame = std::move(frame);
    return showMolecule(m_frame, m_camera);
}

void View::setCamera(const ViewCamera& camera) {
    m_camera = camera;
    if (m_frame) {
        showMolecule(m_frame, m_camera);
    }
}

void View::settle() {
    if (m_frame && !m_dragging && m_camera.quality != m_frame->quality) {
        ViewCamera c = m_camera;
        c.quality = m_frame->quality;
        setCamera(c);
    }
}

bool View::eventFilter(QObject* watched, QEvent* event) {
    if (watched != m_window || !m_frame) {
        return QObject::eventFilter(watched, event);
    }
    switch (event->type()) {
    case QEvent::MouseButtonPress: {
        auto* e = static_cast<QMouseEvent*>(event);
        if (e->button() == Qt::LeftButton || e->button() == Qt::RightButton) {
            m_dragging = true;
            m_last = e->position();
            return true;
        }
        break;
    }
    case QEvent::MouseMove: {
        auto* e = static_cast<QMouseEvent*>(event);
        if (m_dragging) {
            const QPointF d = e->position() - m_last;
            m_last = e->position();
            ViewCamera c = rotatedCamera(m_camera, d.x(), d.y(), m_window->height());
            c.quality = molecule::RenderQuality::Fast;
            setCamera(c);
            return true;
        }
        break;
    }
    case QEvent::MouseButtonRelease:
        if (m_dragging) {
            m_dragging = false;
            settle();
            return true;
        }
        break;
    case QEvent::MouseButtonDblClick:
        m_dragging = false;
        setCamera(frameCamera(*m_frame));
        return true;
    case QEvent::Wheel: {
        auto* e = static_cast<QWheelEvent*>(event);
        const double steps = e->angleDelta().y() / 120.0;
        if (steps != 0.0) {
            ViewCamera c = zoomedCamera(m_camera, steps);
            c.quality = molecule::RenderQuality::Fast;
            setCamera(c);
            m_settle->start();
        }
        return true;
    }
    default:
        break;
    }
    return QObject::eventFilter(watched, event);
}

View::~View() { delete m_window; }  // QPointer: null if a window container deleted it already

void View::countFrame() { emit frameRendered(++m_frames); }

ViewChoice createView(Renderer requested) {
    QString reason = QStringLiteral("requested");
    if (requested != Renderer::Software) {
        if (std::unique_ptr<View> view = createVulkanView(&reason)) {
            return {std::move(view), QString()};
        }
    }
    return {createSoftwareView(reason), reason};
}

std::unique_ptr<View> createSoftwareView(const QString& reason) {
    return std::make_unique<SoftwareView>(reason);
}

}  // namespace mdbridge::render
