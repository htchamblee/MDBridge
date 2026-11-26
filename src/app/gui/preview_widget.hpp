#pragma once
// The figure canvas of the main window: MDATK v0.0.4.1's WYSIWYG layout editor (mdatk/gui/canvas.py)
// on a worker thread. The project's figure is built with the example dataset and watermarked
// (app/commands/project_view.hpp previewFigure) at the widget's size whenever the content or the
// size changes; its structure panel is the 3D view (owner decision D4 as amended).
//
//   drag in the structure panel     rotates the model (orientation::drag, 0.5° per pixel); in map
//   (below its title strip, 12 %)   mode it moves the density map, with Shift it rotates the map
//   drag anywhere else (the title   moves the panel: the slot under the pointer is highlighted
//   strip, Shift held)              with its type label; the release swaps the two panels
//   double-click                    opens the slot's editor (slotActivated)
//   wheel over the structure        zooms (x 1.15 per notch; zoomRequested)
//
// The keyboard (D10; MDATK's canvas has no keyboard use): the canvas takes the focus with Tab or a
// click (a focus frame shows it); then the arrow keys turn the structure as a drag of 10 pixels
// would (5°; with Shift 2 pixels, 1°; keyRotated) and + / - zoom it as a wheel notch (zoomRequested);
// Esc cancels a drag.
//
// A press becomes a drag after 4 pixels (|dx| + |dy|). Hit areas are each panel's drawn extent
// (Figure::tightBox of its axes and the twin axes in its cell), an empty slot's grid cell.
//
// AP-20: the figure is built and drawn on a worker, never on the interface thread. Rotation, zoom
// and map moves change only the structure panel of the figure already built (StructurePanel
// setRotation / setZoom / nudgeDensity) and recompose its moving artists over the static
// background (MDATK's blitting), also on the worker: requests made while it works are queued and
// applied in order, so the picture keeps up and no step is lost. While dragging the panel draws at
// fast quality (and its 3D layer with setFastStageRenderer's renderer: the GPU one when there is
// a GPU); full quality follows on release, and 250 ms after the last slider change (MDATK).

#include <array>
#include <cstdint>
#include <deque>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <QImage>
#include <QRectF>
#include <QTimer>
#include <QWidget>

#include "background.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/rigid_transform.hpp"

namespace mdbridge::render::panels {
class StageRenderer;
}

namespace mdbridge::app {

// What the preview draws: the project, its example dataset and the watermark.
struct PreviewContent {
    core::Project project;
    core::Dataset example{{}, {}};
    std::string watermark;
};

struct CanvasScene;

class PreviewWidget : public QWidget {
    Q_OBJECT

public:
    explicit PreviewWidget(QWidget* parent = nullptr);
    ~PreviewWidget() override;

    // Draws `content` (nullptr: the hint to open a project). The canvas's angles and zoom become
    // the project's (MDATK's set_content).
    void setContent(std::shared_ptr<const PreviewContent> content);
    // A message instead of a figure (e.g. "Cannot draw this layout: ...").
    void setMessage(const QString& message);

    // The figure as last drawn (null before the first one), and why the last drawing failed.
    QImage image() const { return m_image; }
    QString problem() const { return m_problem; }
    // The warnings the last drawing logged ("WARNING: ..."; layout warnings, centre-on misses).
    const std::vector<std::string>& warnings() const { return m_warnings; }
    // The image shows the current content at the current size.
    bool current() const { return m_content && !m_image.isNull() && m_drawn == m_generation; }
    // Nothing more to draw: no content (a message), or the current content drawn or failed, and no
    // interactive change waiting or being drawn (no full-quality redraw due either).
    bool settled() const;

    QSize sizeHint() const override { return {640, 360}; }

    // ---- the editor ----

    // The slots' hit areas in widget pixels (after the first drawing; empty before) and the slot of
    // the structure panel (-1: none that draws).
    std::vector<QRectF> slotRects() const;
    int structureSlot() const { return m_structureSlot; }
    // The slot under a widget point (-1: none).
    int slotAt(const QPointF& point) const;
    // The figure's top-left corner in the widget.
    QPointF figureOffset() const;

    // Map mode (the Structure view page's "Move the map with the mouse").
    void setMapMode(bool on);
    bool mapMode() const { return m_mapMode; }
    // The model orientation from the sliders: the structure panel only, at fast quality, full
    // quality 250 ms after the last change (MDATK's set_rotation(fast=True)). Angles that round to
    // the ones shown (3 decimals, as stored) change nothing, so the canvas's own rotations come back
    // without a step.
    void setRotation(const std::array<double, 3>& angles);
    std::array<double, 3> angles() const { return m_angles; }
    // A new zoom: the structure panel at once (AP-20); the window rebuilds the figure later.
    void setZoom(double zoom);
    // The 3D layer's renderer while dragging (nullptr: the panel's own).
    void setFastStageRenderer(std::shared_ptr<const render::panels::StageRenderer> renderer);

    // The structure panel's model (frame 0 of the example dataset: the reference frame once it is
    // loaded): its centre and its residues in place (StructurePanel::modelPoints, the anchor of a
    // placement made on it); nullopt / empty without a structure panel that draws.
    std::optional<core::Vec3> structureCentre() const { return m_centre; }
    const std::vector<std::array<double, 3>>& structureModel() const { return m_model; }
    // The structure panel shows a density map (map mode drags move it).
    bool structureHasDensity() const { return m_hasDensity; }

    // A drag in progress: "rotate", "map", "map-rotate", "move" ("" without one; a press not yet
    // moved 4 pixels is none).
    QString dragMode() const { return m_dragging && m_press ? m_press->mode : QString(); }
    // The slot a panel drag would drop on (-1: none).
    int dragTarget() const { return m_dragTarget; }
    // Ends a drag without its result (Esc): a panel drag drops nothing; rotations and map moves
    // already made stay with the drag's step, which the window undoes (dragCancelled).
    void cancelDrag();

signals:
    // A drawing finished (successfully or not: problem()).
    void drawn();
    // A drag that changes the project began / ended (rotate, map moves: one undo step each);
    // dragFinished comes once every change of the drag has been reported (mapMoved).
    void dragStarted(const QString& mode);
    void dragFinished();
    void dragCancelled();
    // The model was rotated by a drag (each motion; the angles unrounded).
    void rotated(const std::array<double, 3>& angles);
    // The density map was moved by a drag: its new placement (12 numbers, model -> map).
    void mapMoved(const std::vector<double>& transform);
    // A panel was dropped on another slot.
    void swapRequested(int from, int to);
    // The wheel over the structure panel: multiply the zoom by `factor`.
    void zoomRequested(double factor);
    // A double-click on slot `slot`.
    void slotActivated(int slot);
    // An arrow key turned the model (the angles unrounded; the canvas shows them already).
    void keyRotated(const std::array<double, 3>& angles);

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseDoubleClickEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;
    void keyPressEvent(QKeyEvent* event) override;
    void focusInEvent(QFocusEvent* event) override;
    void focusOutEvent(QFocusEvent* event) override;
    void leaveEvent(QEvent* event) override;

public:
    // One change of the structure panel of the built figure (applied in order on the worker).
    struct Op {
        enum class Kind { Rotate, Zoom, MapShift, MapRotate, Fast, Full };
        Kind kind = Kind::Rotate;
        std::array<double, 3> angles{};
        double zoom = 1.0;
        double dx = 0.0, dy = 0.0;
    };

private:
    struct Press {
        QPointF at;
        int slot = -1;
        QString mode;
    };

    void schedule();
    void start();
    void startInteraction();
    void runNext();
    void queue(const Op& op);
    void updateCursor(const QPointF& at);
    void drawDrag(QPainter& painter);
    void drawFocus(QPainter& painter);
    void settleNow();
    void maybeFinish();

    std::shared_ptr<const PreviewContent> m_content;
    QImage m_image;
    QString m_problem;
    std::vector<std::string> m_warnings;
    QString m_message;
    std::uint64_t m_generation = 0;
    std::uint64_t m_drawn = 0;
    bool m_busy = false;      // a build or an interactive change runs on the worker
    bool m_pending = false;   // a build is due once the worker is free
    QTimer m_debounce;
    QTimer m_settle;

    // the figure built last (owned by the worker while m_busy) and what the interface knows of it
    std::shared_ptr<CanvasScene> m_scene;
    std::deque<Op> m_ops;
    std::vector<QRectF> m_slots;  // figure pixels, y down
    QSizeF m_figureSize;
    int m_structureSlot = -1;
    std::vector<std::string> m_panelTypes;
    std::optional<core::Vec3> m_centre;
    std::vector<std::array<double, 3>> m_model;
    bool m_hasDensity = false;
    std::shared_ptr<const render::panels::StageRenderer> m_fastRenderer;

    std::array<double, 3> m_angles{};
    double m_zoom = 1.0;
    bool m_mapMode = false;
    std::optional<Press> m_press;
    QPointF m_last;
    QPointF m_dragPos;
    bool m_dragging = false;
    bool m_changeDrag = false;  // dragStarted was emitted for the current drag
    bool m_finishPending = false;  // dragFinished waits for the drag's changes to be drawn
    bool m_dropMoves = false;      // a cancelled drag: its map moves still being drawn are not reported
    int m_dragTarget = -1;

    Background m_background{this};
};

}  // namespace mdbridge::app
