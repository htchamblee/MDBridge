#pragma once
// The 3D structure view. Vulkan 1.3 when the machine has it; otherwise the software view, which
// draws fewer representations and effects but takes the same place in the panels. Without a
// molecule (setMoleculeFrame) both draw milestone 1's placeholder scene.

#include <array>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

#include <QImage>
#include <QObject>
#include <QPointF>
#include <QPointer>
#include <QString>
#include <QWindow>

class QTimer;

#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"

namespace mdbridge::render {

// What the 3D view shows: the molecule scene's primitives of one frame (MoleculeScene::frame), the
// view rotation, the data window (Angstrom; the view widens one axis to the window's aspect, so
// atoms stay round), the depth range of the depth cue and the shading, over `background` (RGB
// 0..1). The view draws it at the window's size, as the structure panel's 3D layer would.
struct MoleculeFrame {
    std::vector<molecule::Colored> items;
    std::array<molecule::Vec3, 3> rotation{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    std::array<double, 2> xlim{-1.0, 1.0};
    std::array<double, 2> ylim{-1.0, 1.0};
    std::array<double, 2> zrange{-1.0, 1.0};
    molecule::Shading shading;
    molecule::RenderQuality quality = molecule::RenderQuality::High;  // Fast while dragging
    std::array<float, 3> background{0.102f, 0.118f, 0.153f};
};

// What the mouse changes over a MoleculeFrame: the view rotation (rows), the zoom (> 1 magnifies
// the data window about its centre) and the quality (fast while the user drags).
struct ViewCamera {
    std::array<molecule::Vec3, 3> rotation{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    double zoom = 1.0;
    molecule::RenderQuality quality = molecule::RenderQuality::High;
};

// The camera a frame starts with: its rotation and quality, no zoom.
ViewCamera frameCamera(const MoleculeFrame& frame);
// The viewport that shows `frame` through `camera` in a target of width x height pixels: the data
// window scaled by 1 / zoom about its centre, widened along one axis to the target's aspect.
molecule::Viewport cameraViewport(const MoleculeFrame& frame, const ViewCamera& camera, int width, int height);
// The viewport that shows `frame` in a target of width x height pixels: the data window centred
// and widened along one axis to the target's aspect (cameraViewport with frameCamera(frame)).
molecule::Viewport fittedViewport(const MoleculeFrame& frame, int width, int height);
// The camera after a mouse drag of (dx, dy) pixels (y down) in a view `height` pixels tall: the
// model turns about the screen's vertical and horizontal axes, half a turn per view height, its
// front following the pointer.
ViewCamera rotatedCamera(const ViewCamera& camera, double dx, double dy, double height);
// The camera after `steps` wheel steps (120 units of angleDelta each; positive zooms in): zoom
// times 1.15^steps, kept within [0.05, 50].
ViewCamera zoomedCamera(const ViewCamera& camera, double steps);

enum class Renderer { Auto, Vulkan, Software };

// "auto", "vulkan" or "software"
std::optional<Renderer> rendererFromName(std::string_view name);
const char* rendererName(Renderer renderer);

class View : public QObject {
    Q_OBJECT

public:
    ~View() override;

    // The window that shows the view (embed it with QWidget::createWindowContainer). The view
    // deletes it unless something else did first.
    QWindow* window() const { return m_window; }
    Renderer renderer() const { return m_renderer; }
    // One line for the status bar, e.g. "Vulkan 1.3.275 · AMD Radeon RX 6600".
    QString description() const { return m_description; }
    int frames() const { return m_frames; }

    // The last frame as an image; null if it cannot be read back.
    virtual QImage grabFrame() = 0;

    // Shows `frame` from the next frame on, with the camera reset to frameCamera(frame) (nullptr:
    // the placeholder again). False when this view cannot draw molecules. Mouse on the view's
    // window: drag (left or right button) to rotate, wheel to zoom, double-click to reset.
    bool setMoleculeFrame(std::shared_ptr<const MoleculeFrame> frame);
    std::shared_ptr<const MoleculeFrame> moleculeFrame() const { return m_frame; }
    const ViewCamera& camera() const { return m_camera; }
    // Shows the frame through `camera` (nothing without a frame).
    void setCamera(const ViewCamera& camera);
    // Whether the last frame drawn showed the current molecule at the window's size (false while
    // the software view still renders it, and for the placeholder).
    virtual bool moleculeShown() const { return false; }

signals:
    void frameRendered(int frames);
    // The view cannot run here (e.g. no Vulkan swapchain over a remote display); its owner should
    // replace it with a software view.
    void failed(const QString& reason);

protected:
    View(QWindow* window, Renderer renderer, QString description);
    void countFrame();
    // Draws `frame` (nullptr: the placeholder) through `camera` from the next frame on; false
    // when the view cannot.
    virtual bool showMolecule(std::shared_ptr<const MoleculeFrame> frame, const ViewCamera& camera) = 0;
    bool eventFilter(QObject* watched, QEvent* event) override;

private:
    void settle();

    QPointer<QWindow> m_window;
    Renderer m_renderer;
    QString m_description;
    int m_frames = 0;
    std::shared_ptr<const MoleculeFrame> m_frame;
    ViewCamera m_camera;
    bool m_dragging = false;
    QPointF m_last;
    QTimer* m_settle = nullptr;  // back to the frame's quality after the wheel stops
};

struct ViewChoice {
    std::unique_ptr<View> view;
    QString fallbackReason;  // why the software view was chosen ("" when it is the Vulkan view)
};

// The view for `requested` (needs a QGuiApplication). Auto and Vulkan both fall back to the
// software view when Vulkan 1.3 is not available, so the program always starts.
ViewChoice createView(Renderer requested);
std::unique_ptr<View> createSoftwareView(const QString& reason);

}  // namespace mdbridge::render
