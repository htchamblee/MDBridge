#pragma once

#include <filesystem>
#include <memory>

#include <QMainWindow>

#include "mdbridge/render/view.hpp"

class QLabel;
class QSplitter;

namespace mdbridge::app {

class ProjectPanel;

// The project in a sidebar next to the 3D view.
class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(render::Renderer renderer, QWidget* parent = nullptr);
    ~MainWindow() override;

    // Shows the project; on failure tells the user (interactive) or prints the reason, and returns false.
    bool openProject(const std::filesystem::path& path, bool interactive = true);

    render::View* view() const { return m_view.get(); }
    // Why the 3D view is the software view ("" when it is the Vulkan view).
    QString fallbackReason() const { return m_fallbackReason; }

signals:
    void viewFrameRendered(int frames);

protected:
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;

private:
    void createMenus();
    void setView(render::ViewChoice choice);
    void useSoftwareView(const QString& reason);
    void chooseProject();
    void showGraphics();
    void showAbout();

    ProjectPanel* m_panel = nullptr;
    QSplitter* m_splitter = nullptr;
    QWidget* m_container = nullptr;
    std::unique_ptr<render::View> m_view;
    QString m_fallbackReason;
    QLabel* m_viewLabel = nullptr;
};

}  // namespace mdbridge::app
