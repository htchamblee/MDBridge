#pragma once
// The 3D structure view. Vulkan 1.3 when the machine has it; otherwise the software view, which
// draws fewer representations and effects but takes the same place in the panels. Milestone 1
// draws a placeholder scene in both.

#include <memory>
#include <optional>
#include <string_view>

#include <QImage>
#include <QObject>
#include <QPointer>
#include <QString>
#include <QWindow>

namespace mdbridge::render {

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
