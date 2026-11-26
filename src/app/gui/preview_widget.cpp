#include "preview_widget.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>

#include "log_capture.hpp"
#pragma push_macro("slots")
#undef slots
#include "mdbridge/core/panels.hpp"
#pragma pop_macro("slots")
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure_qt.hpp"
#include "mdbridge/render/orientation.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "project_view.hpp"

namespace mdbridge::app {

namespace {
constexpr int kDebounceMs = 120;      // MDATK's canvas: schedule_render(delay_ms=120)
constexpr int kSettleMs = 250;        // MDATK's _schedule_settle(delay_ms=250)
constexpr double kTitleStrip = 0.12;  // top fraction of a slot that drags the structure panel
constexpr double kDragThreshold = 4;  // pixels before a press becomes a drag
constexpr double kDegreesPerPixel = 0.5;
constexpr double kWheelFactor = 1.15;

// round(x, 3): how visual.rotation is stored
double round3(double x) { return std::round(x * 1000.0) / 1000.0; }

}  // namespace

// The figure built last; used by one thread at a time (the worker while PreviewWidget::m_busy).
struct CanvasScene {
    render::panels::BuiltFigure built;
    const render::FontSet* fonts = nullptr;
    std::size_t pos = 0;
    render::panels::StructurePanel* structure = nullptr;
    QImage background;  // everything but the structure panel's moving artists
    qreal dpr = 1.0;
    std::array<double, 3> angles{};

    // The structure panel's moving artists over the background.
    QImage compose() {
        if (structure == nullptr) {
            return background;
        }
        structure->update(built.figure, pos);
        built.figure.draw(*fonts);
        QImage out = render::composeLayers(background, render::buildDisplayList(built.figure, *fonts, render::DisplayLayer::Animated));
        out.setDevicePixelRatio(dpr);
        return out;
    }
};

namespace {

struct BuildResult {
    std::shared_ptr<CanvasScene> scene;
    QImage image;
    QString problem;
    std::vector<std::string> warnings;
    std::vector<QRectF> rects;  // figure device pixels, y down
    QSizeF size;
    int structureSlot = -1;
    std::vector<std::string> panelTypes;
    std::optional<core::Vec3> centre;
    std::vector<std::array<double, 3>> model;
    bool hasDensity = false;
};

// MDATK's _compute_slots: a panel's drawn extent (its axes and those of its cell), an empty slot's
// grid cell; figure device pixels with y down.
std::vector<QRectF> figureSlotRects(const core::Project& project, const render::panels::BuiltFigure& built) {
    std::vector<QRectF> out;
    const core::LayoutPreset* preset = core::findLayoutPreset(project.layout.preset);
    if (preset == nullptr) {
        return out;
    }
    const double cw = built.figure.widthPx(), ch = built.figure.heightPx();
    std::vector<double> heights(preset->heightRatios.begin(), preset->heightRatios.end());
    if (heights.empty()) heights.assign(static_cast<std::size_t>(preset->rows), 1.0);
    const double total = std::accumulate(heights.begin(), heights.end(), 0.0);
    std::vector<double> edges{0.0};
    for (const double h : heights) edges.push_back(edges.back() + h);
    for (double& e : edges) e /= total;
    for (std::size_t k = 0; k < core::presetSlots(*preset).size(); ++k) {
        const core::LayoutSlot& s = core::presetSlots(*preset)[k];
        const auto r0 = static_cast<std::size_t>(s.row0), r1 = static_cast<std::size_t>(s.row1);
        QRectF rect(cw * s.col0 / preset->cols, ch * edges[r0], cw * (s.col1 - s.col0) / preset->cols,
                    ch * (edges[r1] - edges[r0]));
        if (k < built.panels.size()) {
            const render::panels::Panel& panel = *built.panels[k];
            std::vector<std::size_t> axes{panel.mainAxes()};
            if (panel.twinAxes()) axes.push_back(*panel.twinAxes());
            std::vector<render::Box> boxes;
            for (const std::size_t a : axes) {
                if (const std::optional<render::Box> b = built.figure.tightBox(a); b && b->width() > 0) {
                    boxes.push_back(*b);
                }
            }
            if (!boxes.empty()) {
                const render::Box u = render::unionOf(boxes);
                rect = QRectF(u.x0, ch - u.y1, u.width(), u.height());
            }
        }
        out.push_back(rect);
    }
    return out;
}

}  // namespace

PreviewWidget::PreviewWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(160, 90);
    setMouseTracking(true);
    setFocusPolicy(Qt::ClickFocus);
    m_message = tr("The figure preview appears here once a project is open.");
    m_debounce.setSingleShot(true);
    m_debounce.setInterval(kDebounceMs);
    connect(&m_debounce, &QTimer::timeout, this, &PreviewWidget::start);
    m_settle.setSingleShot(true);
    m_settle.setInterval(kSettleMs);
    connect(&m_settle, &QTimer::timeout, this, &PreviewWidget::settleNow);
}

PreviewWidget::~PreviewWidget() = default;

bool PreviewWidget::settled() const {
    if (!m_content) return true;
    return m_drawn == m_generation && !m_busy && m_ops.empty() && !m_settle.isActive() && !m_finishPending;
}

void PreviewWidget::setContent(std::shared_ptr<const PreviewContent> content) {
    m_content = std::move(content);
    m_message.clear();
    m_problem.clear();
    // rotations and zooms come with the content; map moves not yet drawn (and reported) still apply
    std::erase_if(m_ops, [](const Op& op) { return op.kind != Op::Kind::MapShift && op.kind != Op::Kind::MapRotate; });
    m_settle.stop();
    ++m_generation;
    if (m_content) {
        const std::vector<double>& r = m_content->project.visual.rotation;
        m_angles = r.size() == 3 ? std::array<double, 3>{r[0], r[1], r[2]} : std::array<double, 3>{};
        m_zoom = m_content->project.visual.zoom;
        if (m_image.isNull()) m_message = tr("Building preview ...");
        schedule();
    } else {
        m_image = QImage();
    }
    update();
}

void PreviewWidget::setMessage(const QString& message) {
    m_content.reset();
    m_image = QImage();
    m_ops.clear();
    m_settle.stop();
    m_slots.clear();
    m_structureSlot = -1;
    m_centre.reset();
    m_model.clear();
    m_hasDensity = false;
    ++m_generation;
    m_message = message;
    update();
}

void PreviewWidget::resizeEvent(QResizeEvent* event) {
    QWidget::resizeEvent(event);
    if (m_content) {
        ++m_generation;
        schedule();
    }
}

void PreviewWidget::schedule() { m_debounce.start(); }

void PreviewWidget::setFastStageRenderer(std::shared_ptr<const render::panels::StageRenderer> renderer) {
    m_fastRenderer = std::move(renderer);
}

void PreviewWidget::start() {
    if (!m_content) {
        return;
    }
    if (m_busy) {
        m_pending = true;  // built again when the running job ends
        return;
    }
    m_pending = false;
    const qreal dpr = devicePixelRatioF();
    const int width = static_cast<int>(std::lround(this->width() * dpr));
    const int height = static_cast<int>(std::lround(this->height() * dpr));
    if (width < kPreviewMinimumSize || height < kPreviewMinimumSize) {
        return;  // MDATK draws nothing in a space under 50 pixels
    }
    m_busy = true;
    m_background.run(
        [content = m_content, width, height, dpr, fast = m_fastRenderer]() {
            BuildResult r;
            const ThreadLogCapture capture(&r.warnings);
            try {
                auto scene = std::make_shared<CanvasScene>();
                scene->fonts = &render::bundledFontSet();
                scene->dpr = dpr;
                scene->built = previewFigure(content->project, content->example, *scene->fonts, width, height,
                                             content->watermark);
                const std::size_t frames = content->example.frameCount();
                scene->pos = frames > 0 ? previewPosition(frames) : 0;
                render::panels::BuiltFigure& built = scene->built;
                for (std::size_t k = 0; k < built.panels.size(); ++k) {
                    auto* st = dynamic_cast<render::panels::StructurePanel*>(built.panels[k].get());
                    if (st != nullptr && st->ok()) {
                        scene->structure = st;
                        r.structureSlot = static_cast<int>(k);
                        break;
                    }
                }
                for (const core::PanelSpec& p : content->project.layout.panels) r.panelTypes.push_back(p.type);
                r.rects = figureSlotRects(content->project, built);
                r.size = QSizeF(built.figure.widthPx(), built.figure.heightPx());
                if (render::panels::StructurePanel* st = scene->structure) {
                    st->setFastRenderer(fast);
                    std::vector<render::panels::ArtistRef> moving = st->rotationArtists();
                    for (const render::panels::ArtistRef& a : st->animatedArtists()) moving.push_back(a);
                    render::panels::markAnimated(built.figure, moving);
                    scene->background = render::renderImage(
                        render::buildDisplayList(built.figure, *scene->fonts, render::DisplayLayer::Static));
                    scene->background.setDevicePixelRatio(dpr);
                    r.image = scene->compose();
                    r.centre = core::Vec3{st->center()[0], st->center()[1], st->center()[2]};
                    r.model = st->modelPoints(0);
                    r.hasDensity = st->hasDensity();
                } else {
                    scene->background = render::renderImage(render::buildDisplayList(built.figure, *scene->fonts));
                    scene->background.setDevicePixelRatio(dpr);
                    r.image = scene->background;
                }
                const std::vector<double>& rot = content->project.visual.rotation;
                if (rot.size() == 3) scene->angles = {rot[0], rot[1], rot[2]};
                r.scene = std::move(scene);
            } catch (const std::exception& e) {
                r.problem = tr("Cannot draw this layout: %1").arg(QString::fromUtf8(e.what()));
            }
            return r;
        },
        [this, generation = m_generation](BuildResult r) {
            m_busy = false;
            if (generation == m_generation) {
                m_drawn = generation;
                m_scene = std::move(r.scene);
                m_image = r.image;
                m_problem = r.problem;
                m_message = r.problem;
                m_warnings = std::move(r.warnings);
                const qreal ratio = devicePixelRatioF();
                m_slots.clear();
                for (const QRectF& s : r.rects) {
                    m_slots.emplace_back(s.x() / ratio, s.y() / ratio, s.width() / ratio, s.height() / ratio);
                }
                m_figureSize = r.size / ratio;
                m_structureSlot = m_scene ? r.structureSlot : -1;
                m_panelTypes = std::move(r.panelTypes);
                m_centre = r.centre;
                m_model = std::move(r.model);
                m_hasDensity = r.hasDensity;
                if (m_problem.isEmpty() && m_scene && m_scene->structure != nullptr) {
                    // what changed since this build began (a drag, the wheel) applies to the new figure
                    const std::array<double, 3> built = m_scene->angles;
                    if (built != m_angles) m_ops.push_front(Op{Op::Kind::Rotate, m_angles, 0.0, 0.0, 0.0});
                }
                update();
                emit drawn();
            } else {
                m_pending = true;
            }
            if (m_pending) {
                start();
            } else {
                runNext();
            }
            maybeFinish();
        });
}

// ---- interactive changes ----------------------------------------------------------------------

void PreviewWidget::queue(const Op& op) {
    if (!m_scene || m_scene->structure == nullptr) {
        return;
    }
    // consecutive rotations or zooms: only the last one matters
    if (!m_ops.empty() && m_ops.back().kind == op.kind && (op.kind == Op::Kind::Rotate || op.kind == Op::Kind::Zoom)) {
        m_ops.back() = op;
    } else {
        m_ops.push_back(op);
    }
    runNext();
}

void PreviewWidget::runNext() {
    if (m_busy || m_ops.empty() || !m_scene || m_drawn != m_generation) {
        return;
    }
    m_busy = true;
    std::vector<Op> ops(m_ops.begin(), m_ops.end());
    m_ops.clear();
    struct Result {
        QImage image;
        std::optional<std::vector<double>> moved;  // the last map placement
        std::string problem;
    };
    m_background.run(
        [scene = m_scene, ops = std::move(ops)]() {
            Result r;
            try {
                render::panels::StructurePanel& st = *scene->structure;
                render::Figure& fig = scene->built.figure;
                for (const Op& op : ops) {
                    switch (op.kind) {
                    case Op::Kind::Rotate:
                        scene->angles = op.angles;
                        st.setRotation(fig, op.angles);
                        break;
                    case Op::Kind::Zoom:
                        st.setZoom(fig, op.zoom);
                        break;
                    case Op::Kind::Fast:
                        st.setFast(true);
                        break;
                    case Op::Kind::Full:
                        st.setFast(false);
                        st.setRotation(fig, scene->angles);  // the shadow at full quality
                        break;
                    case Op::Kind::MapShift: {
                        // widget pixels -> Å in the panel's data coordinates (screen y points down)
                        const render::Box bb = fig.axesBox(st.mainAxes());
                        const render::Axes& ax = fig.axes[st.mainAxes()];
                        const double sx = ax.xaxis.limits.hi - ax.xaxis.limits.lo;
                        const double sy = ax.yaxis.limits.hi - ax.yaxis.limits.lo;
                        const double dpr = scene->dpr;
                        if (auto moved = st.nudgeDensity(fig, {op.dx * dpr * sx / std::max(bb.width(), 1.0),
                                                               -op.dy * dpr * sy / std::max(bb.height(), 1.0)})) {
                            r.moved = std::move(moved);
                        }
                        break;
                    }
                    case Op::Kind::MapRotate: {
                        const render::Mat3 rot = render::multiply(render::rotationY(op.dx * kDegreesPerPixel),
                                                                  render::rotationX(op.dy * kDegreesPerPixel));
                        if (auto moved = st.nudgeDensity(fig, {0.0, 0.0}, rot)) {
                            r.moved = std::move(moved);
                        }
                        break;
                    }
                    }
                }
                r.image = scene->compose();
            } catch (const std::exception& e) {
                r.problem = e.what();
            }
            return r;
        },
        [this, generation = m_generation](Result r) {
            m_busy = false;
            if (generation == m_generation) {
                if (!r.image.isNull()) m_image = r.image;
                update();
                if (r.moved && !m_dropMoves) emit mapMoved(*r.moved);
            }
            if (m_pending) {
                start();
            } else {
                runNext();
            }
            maybeFinish();
        });
}

void PreviewWidget::maybeFinish() {
    if (m_finishPending && !m_busy && m_ops.empty()) {
        m_finishPending = false;
        emit dragFinished();
    }
}

void PreviewWidget::setRotation(const std::array<double, 3>& angles) {
    const std::array<double, 3> shown{round3(m_angles[0]), round3(m_angles[1]), round3(m_angles[2])};
    if (angles == shown || angles == m_angles) {
        return;
    }
    m_angles = angles;
    queue(Op{Op::Kind::Fast});
    queue(Op{Op::Kind::Rotate, m_angles});
    m_settle.start();
}

void PreviewWidget::setZoom(double zoom) {
    if (zoom == m_zoom) {
        return;
    }
    m_zoom = zoom;
    Op op{Op::Kind::Zoom};
    op.zoom = zoom;
    queue(op);
}

void PreviewWidget::settleNow() {
    m_settle.stop();
    queue(Op{Op::Kind::Full, m_angles});
}

void PreviewWidget::setMapMode(bool on) {
    m_mapMode = on;
    updateCursor(mapFromGlobal(QCursor::pos()));
}

// ---- geometry ---------------------------------------------------------------------------------

QPointF PreviewWidget::figureOffset() const {
    return {(width() - m_figureSize.width()) / 2.0, (height() - m_figureSize.height()) / 2.0};
}

std::vector<QRectF> PreviewWidget::slotRects() const {
    std::vector<QRectF> out;
    const QPointF o = figureOffset();
    for (const QRectF& s : m_slots) out.push_back(s.translated(o));
    return out;
}

int PreviewWidget::slotAt(const QPointF& point) const {
    if (!current()) {
        return -1;
    }
    const std::vector<QRectF> rects = slotRects();
    for (std::size_t k = 0; k < rects.size(); ++k) {
        const QRectF& s = rects[k];
        if (s.left() <= point.x() && point.x() <= s.right() && s.top() <= point.y() && point.y() <= s.bottom()) {
            return static_cast<int>(k);
        }
    }
    return -1;
}

// ---- the mouse --------------------------------------------------------------------------------

void PreviewWidget::mousePressEvent(QMouseEvent* event) {
    setFocus(Qt::MouseFocusReason);
    const QPointF at = event->position();
    const int slot = slotAt(at);
    if (slot < 0 || event->button() != Qt::LeftButton) {
        QWidget::mousePressEvent(event);
        return;
    }
    const QRectF s = slotRects()[static_cast<std::size_t>(slot)];
    const bool inTitle = at.y() < s.top() + kTitleStrip * s.height();
    const bool shift = (event->modifiers() & Qt::ShiftModifier) != 0;
    const bool inStructure = slot == m_structureSlot && !inTitle;
    QString mode;
    if (inStructure && m_mapMode && m_hasDensity) {
        mode = shift ? QStringLiteral("map-rotate") : QStringLiteral("map");
    } else {
        mode = inStructure && !shift ? QStringLiteral("rotate") : QStringLiteral("move");
    }
    m_press = Press{at, slot, mode};
    m_dropMoves = false;
    m_last = at;
    m_dragging = false;
    m_changeDrag = false;
    m_dragTarget = -1;
    event->accept();
}

void PreviewWidget::mouseMoveEvent(QMouseEvent* event) {
    const QPointF at = event->position();
    if (!m_press) {
        updateCursor(at);
        return;
    }
    if (!m_dragging && std::abs(at.x() - m_press->at.x()) + std::abs(at.y() - m_press->at.y()) < kDragThreshold) {
        return;
    }
    if (!m_dragging) {
        m_dragging = true;
        if (m_press->mode != QStringLiteral("move")) {
            m_changeDrag = true;
            emit dragStarted(m_press->mode);
            queue(Op{Op::Kind::Fast});
        }
        if (m_press->mode == QStringLiteral("rotate")) setCursor(Qt::ClosedHandCursor);
    }
    const double dx = at.x() - m_last.x(), dy = at.y() - m_last.y();
    if (m_press->mode == QStringLiteral("rotate")) {
        m_last = at;
        const render::Vec3d next = render::drag({m_angles[0], m_angles[1], m_angles[2]}, dx, dy, kDegreesPerPixel);
        m_angles = {next[0], next[1], next[2]};
        queue(Op{Op::Kind::Rotate, m_angles});
        emit rotated(m_angles);
    } else if (m_press->mode == QStringLiteral("map") || m_press->mode == QStringLiteral("map-rotate")) {
        m_last = at;
        Op op{m_press->mode == QStringLiteral("map") ? Op::Kind::MapShift : Op::Kind::MapRotate};
        op.dx = dx;
        op.dy = dy;
        queue(op);
    } else {
        m_dragTarget = slotAt(at);
        m_dragPos = at;
        update();
    }
}

void PreviewWidget::mouseReleaseEvent(QMouseEvent* event) {
    if (!m_press) {
        QWidget::mouseReleaseEvent(event);
        return;
    }
    const Press press = *m_press;
    const int target = m_dragTarget;
    const bool dragging = m_dragging;
    const bool changeDrag = m_changeDrag;
    m_press.reset();
    m_dragTarget = -1;
    m_dragging = false;
    m_changeDrag = false;
    if (press.mode == QStringLiteral("move") && dragging && target >= 0 && target != press.slot) {
        emit swapRequested(press.slot, target);
    }
    if (changeDrag) {
        settleNow();  // full-quality shadow and 3D layer after the drag
        m_finishPending = true;
        maybeFinish();
    }
    updateCursor(event->position());
    update();
}

void PreviewWidget::mouseDoubleClickEvent(QMouseEvent* event) {
    const int slot = slotAt(event->position());
    if (slot < 0 || event->button() != Qt::LeftButton) {
        QWidget::mouseDoubleClickEvent(event);
        return;
    }
    m_press.reset();
    m_dragging = false;
    emit slotActivated(slot);
}

void PreviewWidget::wheelEvent(QWheelEvent* event) {
    const int slot = slotAt(event->position());
    if (slot < 0 || slot != m_structureSlot) {
        event->ignore();
        return;
    }
    const double steps = event->angleDelta().y() / 120.0;
    if (steps != 0.0) {
        emit zoomRequested(std::pow(kWheelFactor, steps));
    }
    event->accept();
}

void PreviewWidget::keyPressEvent(QKeyEvent* event) {
    if (event->key() == Qt::Key_Escape && m_press) {
        cancelDrag();
        event->accept();
        return;
    }
    QWidget::keyPressEvent(event);
}

void PreviewWidget::cancelDrag() {
    if (!m_press) return;
    const bool changeDrag = m_changeDrag;
    m_press.reset();
    m_dragTarget = -1;
    m_dragging = false;
    m_changeDrag = false;
    if (changeDrag) {
        m_dropMoves = true;
        std::erase_if(m_ops, [](const Op& op) { return op.kind == Op::Kind::MapShift || op.kind == Op::Kind::MapRotate; });
        emit dragCancelled();
        settleNow();
    }
    unsetCursor();
    update();
}

void PreviewWidget::leaveEvent(QEvent* event) {
    if (!m_press) unsetCursor();
    QWidget::leaveEvent(event);
}

void PreviewWidget::updateCursor(const QPointF& at) {
    const int slot = slotAt(at);
    if (slot < 0) {
        unsetCursor();
        return;
    }
    const QRectF s = slotRects()[static_cast<std::size_t>(slot)];
    const bool rotating = slot == m_structureSlot && at.y() >= s.top() + kTitleStrip * s.height();
    // MDATK's cursors: "grab" (rotate), "all-scroll" (map mode), "move" (panels)
    if (rotating) setCursor(m_mapMode ? Qt::SizeAllCursor : Qt::OpenHandCursor);
    else setCursor(Qt::DragMoveCursor);
}

// ---- drawing ----------------------------------------------------------------------------------

void PreviewWidget::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    if (!m_image.isNull()) {
        painter.fillRect(rect(), m_image.pixelColor(0, 0));  // the figure's page colour around it
        const QSizeF shown = m_image.deviceIndependentSize();
        const QSizeF available = size();
        // a size change shows the last image scaled until the new one is drawn
        const double scale = std::min(available.width() / shown.width(), available.height() / shown.height());
        const QSizeF target = m_drawn == m_generation ? shown : shown * scale;
        const QRectF where(QPointF((available.width() - target.width()) / 2, (available.height() - target.height()) / 2),
                           target);
        painter.setRenderHint(QPainter::SmoothPixmapTransform);
        painter.drawImage(where, m_image);
        if (m_press && m_dragging && m_press->mode == QStringLiteral("move")) {
            drawDrag(painter);
        }
        return;
    }
    painter.fillRect(rect(), palette().window());
    if (!m_message.isEmpty()) {
        painter.setPen(palette().color(QPalette::PlaceholderText));
        painter.drawText(rect().adjusted(12, 12, -12, -12), Qt::AlignCenter | Qt::TextWordWrap, m_message);
    }
}

// MDATK's _draw_drag: the source slot tinted, the target outlined (dashed) and the panel's type at
// the pointer.
void PreviewWidget::drawDrag(QPainter& painter) {
    const std::vector<QRectF> rects = slotRects();
    const auto src = static_cast<std::size_t>(m_press->slot);
    if (src >= rects.size()) return;
    const QColor blue = QColor::fromRgbF(0.16f, 0.47f, 0.84f);
    const auto tint = [&](double alpha) {
        QColor c = blue;
        c.setAlphaF(static_cast<float>(alpha));
        return c;
    };
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rects[src].adjusted(2, 2, -2, -2), tint(0.10));
    if (m_dragTarget >= 0 && static_cast<std::size_t>(m_dragTarget) < rects.size() && m_dragTarget != m_press->slot) {
        const QRectF t = rects[static_cast<std::size_t>(m_dragTarget)].adjusted(3, 3, -3, -3);
        painter.fillRect(t, tint(0.18));
        QPen pen(tint(0.9));
        pen.setWidthF(2.5);
        pen.setDashPattern({8.0 / 2.5, 5.0 / 2.5});  // Qt's dashes are in pen widths
        painter.setPen(pen);
        painter.setBrush(Qt::NoBrush);
        painter.drawRect(t);
    }
    const QString label = src < m_panelTypes.size() ? QString::fromStdString(m_panelTypes[src]) : tr("panel");
    QFont font = painter.font();
    font.setPixelSize(12);
    painter.setFont(font);
    const QFontMetricsF fm(font);
    const QRectF text = fm.tightBoundingRect(label);
    const QPointF at = m_dragPos;
    const QRectF box(at.x() + 10, at.y() + 8, text.width() + 14, text.height() + 10);
    painter.fillRect(box, tint(0.92));
    painter.setPen(Qt::white);
    painter.drawText(QPointF(at.x() + 17, at.y() + 13 + text.height()), label);
}

}  // namespace mdbridge::app
