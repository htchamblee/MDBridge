#pragma once
// The analysis window (milestone 5 W1-A): MDATK v0.0.4.1's layout (mdatk/gui/app.py; owner
// decision D4 as amended) around one core::ProjectDocument and its forms::FormContext.
//
//   header      Open (+ drop-down: recent projects, Demo project…, Browse…), Save, Save as…,
//               Workflow (planner, wizard, PDB clean-up, builder: shells until milestones 6 / 7,
//               D2), Undo / Redo (D5), the main menu (New project (defaults), Close project, Create
//               demo system…, Project log…, Export citations (BibTeX)…, Citations…, Graphics…,
//               About MDBridge)
//   page bar    six page buttons (Ctrl+1 .. Ctrl+6), the canvas hint, Refresh (F5)
//   side bar    520 px: the page title and the page (pages/page_registry.hpp) in a scroll area
//   canvas      the figure preview (the layout editor and its structure panel as the 3D view come
//               with W2-A; the separate 3D view of M1 / M4 is gone, D4)
//   log         jobs::LogView, folded with Ctrl+L
//   action bar  Render (Ctrl+R), Cancel, progress, Open output folder, Open video (the jobs and
//               desktop actions: actions/window_actions.hpp, W2-B)
//   status bar  loading messages and the reference frame's summary
//
// Unsaved changes (AP-13, D6): the title carries the modified marker; New, Open, Close project and
// closing the window ask Save / Discard / Cancel; a recovery copy of a modified project is written
// every 2 minutes (core::writeRecovery) and offered at the next start (offerRecoveries()). A file
// another program changed is noticed (a bar offers to reload it) and Save asks before replacing
// it. Questions, messages and file dialogs go through prompter() and forms::runFileDialog, so the
// GUI tests answer them.

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <QMainWindow>
#include <QPointer>
#include <QTimer>

#include "background.hpp"
#include "gui/dialogs/text_dialog.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/pages/page_registry.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/view.hpp"

class QAbstractButton;
class QAction;
class QButtonGroup;
class QFileSystemWatcher;
class QFrame;
class QLabel;
class QMenu;
class QProgressBar;
class QPushButton;
class QStackedWidget;
class QToolButton;

namespace mdbridge::app {

namespace jobs {
class JobRunner;
class LogView;
}  // namespace jobs

namespace problems {
class ProblemsPanel;
}

class PreviewWidget;
class Prompter;
class TextDialog;
class WindowActions;

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    // `renderer`: the 3D renderer the structure panel will draw with (W2-A); kept for it and shown
    // by Graphics….
    explicit MainWindow(render::Renderer renderer = render::Renderer::Auto, QWidget* parent = nullptr);
    ~MainWindow() override;

    // ---- start-up (`mdbridge gui [PROJECT]`) ----

    // Opens the project named on the command line. When it does not load, the window keeps the
    // default settings, says so on stderr (MDATK's `warning: could not load <path>: <reason>`) and,
    // as MDATK 0.1.0 does (AP-12), in its log (unfolded) and status bar. A parameter log opens as an
    // unsaved project. Returns false when the file did not load.
    bool startWith(const std::filesystem::path& path);
    // Offers the recovery copies left by a session that ended without saving (D6), newest first.
    void offerRecoveries();

    // ---- the user's actions (each asks what the menu item asks) ----

    void newProject();                                 // Ctrl+N
    void closeProject();                               // Ctrl+W
    void chooseAndOpen();                              // Ctrl+O, Browse…
    // Opens `path` (recent list, drag and drop, demo); asks about unsaved changes first when
    // `askUnsaved`. False when it was not opened (refused, cancelled, or it did not load).
    bool openPath(const std::filesystem::path& path, bool askUnsaved = true);
    bool save();                                       // Ctrl+S (asks for a file without one)
    bool saveAs();                                     // Ctrl+Shift+S
    void showProjectLog();
    void exportCitations();
    void showCitations();
    void showAbout();
    void showGraphics();
    void createDemo();
    void refreshPreview();                             // F5
    // Reloads the project file from disk (the external-change bar's Reload).
    void reloadFromDisk();

    // ---- state ----

    core::ProjectDocument& document() { return m_document; }
    const core::ProjectDocument& document() const { return m_document; }
    forms::FormContext& formContext() { return *m_context; }
    jobs::LogView* logView() const { return m_log; }
    jobs::JobRunner* jobRunner() const { return m_jobs; }
    PreviewWidget* preview() const { return m_preview; }
    // The window's actions by objectName ("action.save", "action.open", "action.page.input" ...).
    QAction* action(const QString& name) const;
    // The page buttons' QAbstractButton ("page:<id>" buttons) for W1-D's badges.
    QAbstractButton* pageButton(const QString& id) const;
    // The side bar's pages (pages/page_registry.hpp), in MDATK's order.
    struct Page {
        QString id;
        QString title;
        QWidget* widget = nullptr;
        bool placeholder = false;  // a PlaceholderPage until the page's own is delivered
    };
    const std::vector<Page>& pages() const { return m_pages; }
    // The problems list (W1-D), shown by the header's problems button.
    problems::ProblemsPanel* problemsPanel() const { return m_problems; }

    // The page shown ("input", "groups", ...), and showing one.
    QString currentPage() const;
    bool showPage(const QString& id);
    // Shows `page` and scrolls to its section `section` (objectName "section:<page>.<section>"; ""
    // for the top), as MDATK's open_section. False for an unknown page or section.
    bool openSection(const QString& page, const QString& section = QString());
    // Shows the page of field path `path` (core::problemPage) and the field (W1-D's problems list).
    bool revealField(const QString& path);

    // The questions and messages (tests replace them; the default shows dialogs). Not owned.
    Prompter& prompter() const { return *m_prompter; }
    void setPrompter(Prompter* prompter);

    // An informational line in the log (as MDATK's append_log).
    void appendLog(const QString& text);
    // MDATK's error() / info(): the log line, then the message (prompter()).
    void showError(const QString& text, const QString& detail = QString());
    void showInfo(const QString& text, const QString& detail = QString());

    // The located form errors (typed text a field could not take), one per line in page order; ""
    // when none (Save refuses then: "Some settings are invalid").
    QString formErrorText();

    // The file shown in the title and questions ("" for a new project or a parameter log, which
    // Save does not write).
    std::filesystem::path projectFile() const;

    // The structure and preview data of the project are being loaded.
    bool loading() const { return m_loading; }
    // After a load: the reference frame's structure (for the 3D renderer check and W2-A), and why
    // there is none.
    bool structureLoaded() const { return m_molecule != nullptr; }
    std::shared_ptr<const render::MoleculeFrame> moleculeFrame() const { return m_molecule; }
    // The reference frame as the figure's example data shows it (null until it is loaded).
    std::shared_ptr<const render::StructurePreview> structurePreview() const { return m_structurePreview; }
    // Loads the reference frame again (after cleaning wrote new files).
    void reloadReferenceFrame() { startLoad(); }
    // Render, Render still, Export PDBs, Re-render, Inspect system, Run cleaning, Open output
    // folder / video, Centre on model (W2-B).
    WindowActions* windowActions() const { return m_windowActions; }
    QString loadProblem() const { return m_loadProblem; }
    QString structureSummary() const { return m_structureSummary; }
    void cancelLoad();
    // The preview shows the current settings (no redraw pending or running).
    bool previewSettled() const;
    render::Renderer renderer() const { return m_renderer; }

    // The recovery copy written for the current document ("" when none is written).
    const std::filesystem::path& recoveryFile() const { return m_recoveryFile; }
    // Writes the recovery copy now when the project is modified (the timer's work).
    void writeRecoveryNow();
    // The text dialog opened last (Project log, Citations, About; null once closed).
    TextDialog* lastDialog() const { return m_dialog; }
    // The external-change bar is shown.
    bool externalChangeShown() const;
    // Checks whether the file changed on disk (the file watcher and window activation call it).
    void checkDisk();

signals:
    // The structure and preview data of a project finished loading (or failed to).
    void projectLoaded();
    // The page shown changed.
    void pageChanged(const QString& id);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;
    void closeEvent(QCloseEvent* event) override;
    void changeEvent(QEvent* event) override;
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;

private:
    QAction* makeAction(const QString& name, const QString& text, const QKeySequence& key, const QString& tip);
    void createActions();
    QWidget* createHeader();
    QWidget* createPageBar();
    QWidget* createSideBar();
    QWidget* createActionBar();
    QWidget* createExternalBar();
    void fillOpenMenu();
    void updateTitle();
    void updateUndo();
    void onDocumentChanged(const core::DocumentChange& change);
    // Save / Discard / Cancel when modified (true: go on). `doing`: what is about to happen.
    bool confirmUnsaved(const QString& doing);
    bool saveTo(const std::filesystem::path& target, bool asNew);
    std::filesystem::path askSavePath();
    void afterReplace(const std::filesystem::path& opened);
    void dropRecovery();
    void watchFile();
    void startLoad();
    void schedulePreview();
    void updatePreview();
    QString ask(const QString& kind, const QString& path, const QString& title, const QString& text,
                const QString& detail, const std::vector<std::pair<QString, QString>>& choices,
                const QString& defaultId = QString());
    TextDialog* openTextDialog(const QString& title, const QString& heading, const QString& text, bool monospace);
    struct CitationsOf;

    render::Renderer m_renderer;
    core::ProjectDocument m_document;
    std::unique_ptr<forms::FormContext> m_context;
    Prompter* m_prompter = nullptr;

    std::map<QString, QAction*> m_actions;
    QMenu* m_openMenu = nullptr;
    QButtonGroup* m_pageButtons = nullptr;
    QLabel* m_sideTitle = nullptr;
    QStackedWidget* m_stack = nullptr;
    std::vector<Page> m_pages;
    PreviewWidget* m_preview = nullptr;
    problems::ProblemsPanel* m_problems = nullptr;
    WindowActions* m_windowActions = nullptr;
    jobs::LogView* m_log = nullptr;
    jobs::JobRunner* m_jobs = nullptr;
    QProgressBar* m_progress = nullptr;
    QFrame* m_externalBar = nullptr;
    QLabel* m_externalText = nullptr;
    QFileSystemWatcher* m_watcher = nullptr;

    QPointer<TextDialog> m_dialog;
    QTimer m_recoveryTimer;
    std::filesystem::path m_recoveryFile;
    bool m_externalNoticed = false;
    bool m_closing = false;
    bool m_closeAccepted = false;

    // the reference frame and the preview
    std::shared_ptr<const render::MoleculeFrame> m_molecule;
    std::shared_ptr<const render::StructurePreview> m_structurePreview;
    bool m_loading = false;
    QString m_loadProblem;
    QString m_structureSummary;
    std::uint64_t m_loadGeneration = 0;
    std::shared_ptr<std::atomic<bool>> m_cancel;
    QTimer m_previewTimer;
    QTimer m_reloadTimer;
    std::vector<std::string> m_lastWarnings;  // of the last drawing of the preview
    bool m_warningsKnown = false;
    Background m_background{this};
    core::LogSink m_previousSink;
};

}  // namespace mdbridge::app
