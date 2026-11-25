#pragma once

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <memory>

#include <QMainWindow>

#include "background.hpp"
#include "mdbridge/render/view.hpp"

class QLabel;
class QProgressBar;
class QSplitter;

namespace mdbridge::app {

class PreviewWidget;
class ProjectPanel;

// The project in a sidebar next to the 3D view (the project's reference-frame structure) and the
// read-only figure preview below it. Opening a project shows its summary at once; the structure
// and the preview's example data are loaded on a worker thread, with progress and problems in the
// status bar.
class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(render::Renderer renderer, QWidget* parent = nullptr);
    ~MainWindow() override;

    // Shows the project; on failure tells the user (interactive) or prints the reason, and returns false.
    bool openProject(const std::filesystem::path& path, bool interactive = true);

    render::View* view() const { return m_view.get(); }
    PreviewWidget* preview() const { return m_preview; }
    // Why the 3D view is the software view ("" when it is the Vulkan view).
    QString fallbackReason() const { return m_fallbackReason; }

    // The structure and preview data of the opened project are being loaded.
    bool loading() const { return m_loading; }
    // After a load: the 3D view has the project's structure; else why not ("" while none was tried).
    bool structureLoaded() const { return m_molecule != nullptr; }
    QString loadProblem() const { return m_loadProblem; }
    // e.g. "Reference frame 3: 24 residues, 517 atoms"
    QString structureSummary() const { return m_structureSummary; }
    // Stops a load in progress (its results are dropped).
    void cancelLoad();

signals:
    void viewFrameRendered(int frames);
    // The structure and preview data of a project finished loading (or failed to).
    void projectLoaded();

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
    void startLoad(const std::filesystem::path& path);

    ProjectPanel* m_panel = nullptr;
    QSplitter* m_splitter = nullptr;
    QSplitter* m_right = nullptr;  // the 3D view above the preview
    PreviewWidget* m_preview = nullptr;
    QWidget* m_container = nullptr;
    std::unique_ptr<render::View> m_view;
    QString m_fallbackReason;
    QLabel* m_viewLabel = nullptr;
    QProgressBar* m_progress = nullptr;

    std::shared_ptr<const render::MoleculeFrame> m_molecule;  // kept for a view chosen later
    bool m_loading = false;
    QString m_loadProblem;
    QString m_structureSummary;
    std::uint64_t m_loadGeneration = 0;
    std::shared_ptr<std::atomic<bool>> m_cancel;
    Background m_background{this};
};

}  // namespace mdbridge::app
