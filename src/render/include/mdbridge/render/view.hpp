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
#include <QPointer>
#include <QString>
#include <QWindow>

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

// The viewport that shows `frame` in a target of width x height pixels: the data window centred
// and widened along one axis to the target's aspect.
molecule::Viewport fittedViewport(const MoleculeFrame& frame, int width, int height);

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

    // Shows `frame` from the next frame on (nullptr: the placeholder again). False when this view
    // cannot draw molecules (the software view of milestone 1).
    virtual bool setMoleculeFrame(std::shared_ptr<const MoleculeFrame> frame) {
        (void)frame;
        return false;
    }

signals:
    void frameRendered(int frames);
    // The view cannot run here (e.g. no Vulkan swapchain over a remote display); its owner should
    // replace it with a software view.
    void failed(const QString& reason);

protected:
    View(QWindow* window, Renderer renderer, QString description);
    void countFrame();

private:
    QPointer<QWindow> m_window;
    Renderer m_renderer;
    QString m_description;
    int m_frames = 0;
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
