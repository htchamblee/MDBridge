// mdbridge_gui_window_tests: the analysis window's shell (milestone 5 W1-A): its layout and
// shortcuts (D4, D11), every header and menu action, the questions about unsaved changes (AP-13),
// Save as (AP-11), the start with a file that does not load (AP-12), closing during a job (AP-14),
// recovery copies (D6), the external-change notice, drag and drop, and screenshots of the window
// (one per page, CI artifacts, never compared) and that every page fits the 520 px side bar (also at
// 150 % and 200 %), paths shown with the platform's separators (WN-12), keyboard and screen-reader
// access (D10, WN-13) and where the settings are (D16, WN-14). Qt Test on the offscreen platform; questions,
// messages and file dialogs answered by the harness's RecordingPrompter and a file dialog function.

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include <QAbstractButton>
#include <QAbstractItemView>
#include <QAbstractSlider>
#include <QAbstractSpinBox>
#include <QAccessible>
#include <QComboBox>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QTextEdit>
#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QFontInfo>
#include <QHeaderView>
#include <QLabel>
#include <QLocale>
#include <QMenu>
#include <QMimeData>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QStatusBar>
#include <QScrollBar>
#include <QTableView>
#include <QToolButton>
#include <QTemporaryDir>
#include <QDialog>
#include <QTest>
#include <QUrl>

#include "citations.hpp"
#include "gui/dialogs/text_dialog.hpp"
#include "gui/display_path.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "gui/main_window.hpp"
#include "gui/preview_widget.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/settings.hpp"
#include "mdbridge/version.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;

namespace {

const fs::path kData = core::pathFromUtf8(MDBRIDGE_TEST_DATA);

QString qs(const fs::path& p) { return QString::fromStdString(core::pathToUtf8(p)); }
// A path as the window shows it in its title, recent list, change notice, log, dialogs and status
// bar: the document's absolute path, normalised, with the platform's separators throughout
// (gui::displayPath, WN-12: '\\' on Windows, where the fixture's folder is Qt's '/' path joined by
// std::filesystem's '\\'). The same as qs(p) on Linux and macOS.
QString shown(const fs::path& p) { return app::gui::displayPath(p.lexically_normal()); }
std::string shownText(const fs::path& p) { return shown(p).toStdString(); }

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// A window with the harness's prompter and file dialogs that answer from `files`.
struct Fixture {
    QTemporaryDir tmp;
    fs::path dir;
    guitest::RecordingPrompter prompter;
    QStringList nextFile;  // the next file dialog's answer (empty: cancelled)
    std::vector<app::forms::FileDialogRequest> requests;
    app::forms::FileDialogFunction previous;
    std::optional<app::MainWindow> window;

    Fixture() {
        dir = core::pathFromUtf8(tmp.path().toStdString());
        previous = app::forms::setFileDialogFunction([this](QWidget*, const app::forms::FileDialogRequest& r) {
            requests.push_back(r);
            QStringList reply = nextFile;
            nextFile.clear();
            return reply;
        });
        window.emplace(render::Renderer::Software);
        window->setPrompter(&prompter);
        window->resize(1320, 920);
        window->show();
        if (!QTest::qWaitForWindowExposed(&*window)) qWarning("the window was not exposed");
    }
    ~Fixture() {
        if (window) {
            window->cancelLoad();
            window->setPrompter(nullptr);
            window.reset();
        }
        app::forms::setFileDialogFunction(previous);
    }
    fs::path copy(const char* project, const char* name) {
        const fs::path target = dir / name;
        fs::copy_file(kData / "projects" / project, target, fs::copy_options::overwrite_existing);
        return target;
    }
    std::vector<std::string> log() {
        window->logView()->flush();
        std::vector<std::string> out;
        for (const QString& line : window->logView()->text().split(QLatin1Char('\n'), Qt::SkipEmptyParts)) {
            out.push_back(line.toStdString());
        }
        return out;
    }
    std::string lastLog() {
        const std::vector<std::string> lines = log();
        return lines.empty() ? std::string() : lines.back();
    }
    void answer(const char* kind, const char* choice) {
        prompter.answers = core::Value(core::Table{{kind, std::string(choice)}});
    }
    std::string lastQuestion() const {
        if (prompter.questions.empty()) return {};
        const core::Value& q = prompter.questions.back();
        return q.find("kind")->stringOr() + ":" + q.find("answer")->stringOr();
    }
};

}  // namespace

class WindowTests : public QObject {
    Q_OBJECT

private slots:
    void initTestCase();
    void layoutAndShortcuts();
    void pagesAndSections();
    void titleMarksChanges();
    void newProjectAsks();
    void closeProject();
    void openAndRecentList();
    void openFailure();
    void saveAndSaveAs();
    void saveAsAppendsTheSuffixBeforeAsking();
    void formErrorsRefuseSaving();
    void undoRedo();
    void projectLog();
    void citations();
    void aboutAndGraphics();
    void workflowShells();
    void demoProject();
    void startWithBrokenFile();
    void closingAsks();
    void closingDuringAJob();
    void recoveryCopies();
    void externalChange();
    void dragAndDrop();
    void pagesFitTheSideBar();
    void pathsShownNatively();
    void everyFocusableWidgetIsNamed();
    void everyControlTakesTheKeyboardFocus();
    void tabOrderVisitsEveryField();
    void escapeClosesDialogs();
    void settingsLocations();
    void screenshots();

private:
    QTemporaryDir m_config;
};

void WindowTests::initTestCase() {
    QVERIFY(m_config.isValid());
    // the recent list and recovery copies of these tests stay in a folder of their own (D16)
    qputenv("XDG_CONFIG_HOME", m_config.path().toUtf8());
    qputenv("APPDATA", m_config.path().toUtf8());
}

void WindowTests::layoutAndShortcuts() {
    Fixture f;
    app::MainWindow& w = *f.window;
    // MDATK's layout (D4): header, page bar, side bar of 520 px, preview right, log, action bar
    for (const char* name : {"header", "pagebar", "sidebar", "preview", "log", "actionbar", "progress"}) {
        QVERIFY2(w.findChild<QWidget*>(QString::fromLatin1(name)), name);
    }
    QCOMPARE(w.findChild<QWidget*>(QStringLiteral("sidebar"))->width(), 520);
    QCOMPARE(w.findChild<QProgressBar*>(QStringLiteral("progress"))->format(), QStringLiteral("Idle"));
    QVERIFY(!w.findChild<QWidget*>(QStringLiteral("view3d")));  // no separate 3D view any more (D4)
    QCOMPARE(w.pages().size(), std::size_t{6});
    // D11
    const std::vector<std::pair<const char*, QKeySequence>> keys{
        {"open", QKeySequence::Open},
        {"save", QKeySequence(Qt::CTRL | Qt::Key_S)},
        {"save_as", QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_S)},
        {"render", QKeySequence(Qt::CTRL | Qt::Key_R)},
        {"new", QKeySequence(Qt::CTRL | Qt::Key_N)},
        {"close_project", QKeySequence(Qt::CTRL | Qt::Key_W)},
        {"quit", QKeySequence(Qt::CTRL | Qt::Key_Q)},
        {"undo", QKeySequence(Qt::CTRL | Qt::Key_Z)},
        {"redo", QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_Z)},
        {"refresh", QKeySequence(Qt::Key_F5)},
        {"planner", QKeySequence(Qt::CTRL | Qt::Key_P)},
        {"wizard", QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_W)},
        {"page.input", QKeySequence(Qt::CTRL | Qt::Key_1)},
        {"page.outputs", QKeySequence(Qt::CTRL | Qt::Key_6)},
    };
    for (const auto& [name, key] : keys) {
        QAction* a = w.action(QString::fromLatin1(name));
        QVERIFY2(a, name);
        QCOMPARE(a->shortcut(), key);
    }
    // Ctrl+3 shows the Analyses page, Ctrl+L folds the log
    QTest::keyClick(&w, Qt::Key_3, Qt::ControlModifier);
    QCOMPARE(w.currentPage(), QStringLiteral("analyses"));
    QVERIFY(w.pageButton(QStringLiteral("analyses"))->isChecked());
    QVERIFY(!w.logView()->expanded());  // folded at first, as MDATK's expander
    QTest::keyClick(&w, Qt::Key_L, Qt::ControlModifier);
    QVERIFY(w.logView()->expanded());
    QTest::keyClick(&w, Qt::Key_L, Qt::ControlModifier);
    QVERIFY(!w.logView()->expanded());
    // the action bar's jobs (W2-B, tests/gui/test_actions.cpp): Render is ready, Cancel only while a job runs
    QVERIFY(w.action(QStringLiteral("render"))->isEnabled());
    QVERIFY(!w.action(QStringLiteral("cancel"))->isEnabled());
    // the page buttons carry W1-D's badges
    for (const char* page : {"input", "groups", "analyses", "structure", "layout", "outputs"}) {
        QVERIFY2(w.findChild<QWidget*>(QStringLiteral("badge:") + QString::fromLatin1(page)), page);
    }
    QVERIFY(w.problemsPanel());
}

void WindowTests::pagesAndSections() {
    Fixture f;
    app::MainWindow& w = *f.window;
    for (const app::MainWindow::Page& page : w.pages()) {
        QVERIFY(w.showPage(page.id));
        QCOMPARE(w.currentPage(), page.id);
        QCOMPARE(w.findChild<QLabel*>(QStringLiteral("side.title"))->text(), QStringLiteral("<b>%1</b>").arg(page.title.toHtmlEscaped()));
    }
    QVERIFY(w.openSection(QStringLiteral("analyses"), QStringLiteral("contacts")));
    QCOMPARE(w.currentPage(), QStringLiteral("analyses"));
    QVERIFY(!w.openSection(QStringLiteral("analyses"), QStringLiteral("nope")));
    QVERIFY(!w.openSection(QStringLiteral("nope")));
    QVERIFY(w.revealField(QStringLiteral("output.video.fps")));
    QCOMPARE(w.currentPage(), QStringLiteral("outputs"));
    QVERIFY(w.revealField(QStringLiteral("analyses.contacts.cutoff")));
    QCOMPARE(w.currentPage(), QStringLiteral("analyses"));
    QTRY_COMPARE(QApplication::focusWidget() ? QApplication::focusWidget()->objectName() : QString(),
                 QStringLiteral("analyses.contacts.cutoff"));
    QVERIFY(w.revealField(QStringLiteral("groups.1.name")));  // a list editor's cell
    QCOMPARE(w.currentPage(), QStringLiteral("groups"));
}

void WindowTests::titleMarksChanges() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QCOMPARE(w.windowTitle(), QStringLiteral("unsaved project[*] - MDBridge"));
    QVERIFY(!w.isWindowModified());
    const fs::path p = f.copy("default.mdatk.toml", "p.mdatk.toml");
    QVERIFY(w.openPath(p));
    QCOMPARE(w.windowTitle(), QStringLiteral("%1[*] - MDBridge").arg(shown(p)));
    w.document().set("analyses.contacts.cutoff", 6.0);
    QVERIFY(w.isWindowModified());
    w.document().undo();
    QVERIFY(!w.isWindowModified());
}

void WindowTests::newProjectAsks() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path p = f.copy("default.mdatk.toml", "p.mdatk.toml");
    QVERIFY(w.openPath(p));
    w.document().set("notes", std::string("edited"));
    f.answer("unsaved", "cancel");
    w.action(QStringLiteral("new"))->trigger();
    QCOMPARE(f.lastQuestion(), std::string("unsaved:cancel"));
    QCOMPARE(w.document().path(), p);
    f.answer("unsaved", "discard");
    w.action(QStringLiteral("new"))->trigger();
    QVERIFY(!w.document().hasPath());
    QVERIFY(!w.document().isModified());
    QCOMPARE(f.lastLog(), std::string("New project with default settings (unsaved)."));
    QVERIFY(readFile(p).find("edited") == std::string::npos);  // discarded, the file untouched
    // unmodified: no question
    const std::size_t asked = f.prompter.questions.size();
    w.action(QStringLiteral("new"))->trigger();
    QCOMPARE(f.prompter.questions.size(), asked);
}

void WindowTests::closeProject() {
    Fixture f;
    app::MainWindow& w = *f.window;
    w.action(QStringLiteral("close_project"))->trigger();
    QCOMPARE(f.lastLog(), std::string("Closed the project."));
    QVERIFY(w.openPath(f.copy("default.mdatk.toml", "p.mdatk.toml")));
    w.document().set("notes", std::string("to keep"));
    f.answer("unsaved", "save");
    w.action(QStringLiteral("close_project"))->trigger();
    QCOMPARE(f.lastLog(), std::string("Closed p.mdatk.toml: default settings (unsaved)."));
    QVERIFY(!w.document().hasPath());
    QVERIFY(readFile(f.dir / "p.mdatk.toml").find("to keep") != std::string::npos);  // saved first
}

void WindowTests::openAndRecentList() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path a = f.copy("default.mdatk.toml", "a.mdatk.toml");
    const fs::path b = f.copy("full.mdatk.toml", "b.mdatk.toml");
    f.nextFile = {qs(a)};
    w.action(QStringLiteral("open"))->trigger();
    QCOMPARE(w.document().path(), a);
    QCOMPARE(f.requests.back().title, QStringLiteral("Open project"));
    QCOMPARE(f.lastLog(), "Opened project " + shownText(a));
    f.nextFile = {qs(b)};
    w.action(QStringLiteral("browse"))->trigger();
    QCOMPARE(w.document().path(), b);
    // the drop-down: recent projects (shared settings file, D16), Demo project…, Browse…
    auto* menu = w.findChild<QMenu*>(QStringLiteral("menu.open"));
    QVERIFY(menu);
    emit menu->aboutToShow();
    QStringList texts;
    for (QAction* item : menu->actions()) {
        if (!item->isSeparator()) texts.push_back(item->text());
    }
    QCOMPARE(texts.mid(0, 3), (QStringList{QStringLiteral("Recent projects"), QStringLiteral("b.mdatk.toml"),
                                           QStringLiteral("a.mdatk.toml")}));
    QCOMPARE(texts.mid(texts.size() - 2), (QStringList{QStringLiteral("Demo project…"), QStringLiteral("Browse…")}));
    QAction* recent = menu->findChild<QAction*>(QStringLiteral("recent:1"));
    QVERIFY(recent);
    QCOMPARE(recent->toolTip(), shown(a));
    recent->trigger();
    QCOMPARE(w.document().path(), a);
    // a parameter log opens as an unsaved project (MDATK)
    const fs::path log = f.copy("MDATK_params_log.toml", "MDATK_params_log.toml");
    QVERIFY(w.openPath(log));
    QVERIFY(w.document().isParamsLog());
    QVERIFY(w.projectFile().empty());
    QVERIFY(w.windowTitle().contains(QStringLiteral("settings of a run")));
}

void WindowTests::openFailure() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path broken = f.copy("broken.mdatk.toml", "broken.mdatk.toml");
    QVERIFY(!w.openPath(broken));
    QCOMPARE(f.prompter.dialogs.size(), std::size_t{1});
    QCOMPARE(f.prompter.dialogs[0].find("text")->stringOr(), "Could not open " + shownText(broken));
    QVERIFY(f.lastLog().rfind("ERROR: Could not open " + shownText(broken) + ": ", 0) == 0);
    QVERIFY(!w.document().hasPath());
}

void WindowTests::saveAndSaveAs() {
    Fixture f;
    app::MainWindow& w = *f.window;
    // a new project: Save asks for a file, with MDATK's suggested name
    f.nextFile = {qs(f.dir / "first")};
    w.action(QStringLiteral("save"))->trigger();
    QVERIFY(fs::exists(f.dir / "first.mdatk.toml"));
    QCOMPARE(f.requests.back().title, QStringLiteral("Save project"));
    QVERIFY(f.requests.back().start.endsWith(QStringLiteral("MDATK.mdatk.toml")));
    QVERIFY(!f.requests.back().confirmOverwrite);
    QCOMPARE(w.document().path(), f.dir / "first.mdatk.toml");
    QCOMPARE(f.lastLog(), "Saved project " + shownText(f.dir / "first.mdatk.toml"));
    QVERIFY(!w.isWindowModified());
    // Save writes its own file without a dialog
    w.document().set("notes", std::string("second"));
    const std::size_t dialogs = f.requests.size();
    w.action(QStringLiteral("save"))->trigger();
    QCOMPARE(f.requests.size(), dialogs);
    QVERIFY(readFile(f.dir / "first.mdatk.toml").find("second") != std::string::npos);
    // Save as: a new file, now the window's
    f.nextFile = {qs(f.dir / "copy.mdatk.toml")};
    w.action(QStringLiteral("save_as"))->trigger();
    QCOMPARE(w.document().path(), f.dir / "copy.mdatk.toml");
    // cancelled dialog: nothing happens
    w.action(QStringLiteral("save_as"))->trigger();
    QCOMPARE(w.document().path(), f.dir / "copy.mdatk.toml");
    // the recent list has both, newest first
    const std::vector<fs::path> recent = core::recentProjects();
    QVERIFY(recent.size() >= 2);
    QCOMPARE(recent[0], f.dir / "copy.mdatk.toml");
    QCOMPARE(recent[1], f.dir / "first.mdatk.toml");
}

// AP-11: the question is about the file that is written (the suffix added first)
void WindowTests::saveAsAppendsTheSuffixBeforeAsking() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path existing = f.copy("full.mdatk.toml", "existing.mdatk.toml");
    const std::string before = readFile(existing);
    f.answer("overwrite", "cancel");
    f.nextFile = {qs(f.dir / "existing")};
    QVERIFY(!w.saveAs());
    QCOMPARE(f.lastQuestion(), std::string("overwrite:cancel"));
    QCOMPARE(f.prompter.questions.back().find("path")->stringOr(), shownText(existing));
    QCOMPARE(readFile(existing), before);
    f.answer("overwrite", "replace");
    f.nextFile = {qs(f.dir / "existing")};
    QVERIFY(w.saveAs());
    QVERIFY(readFile(existing) != before);
    QVERIFY(!fs::exists(f.dir / "existing"));
}

void WindowTests::formErrorsRefuseSaving() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("default.mdatk.toml", "p.mdatk.toml")));
    guitest::typeKey(&w, QStringLiteral("inputs.stop"), QStringLiteral("x"));
    QVERIFY(!w.save());
    QCOMPARE(f.prompter.dialogs.back().find("text")->stringOr(), std::string("Some settings are invalid"));
    QCOMPARE(f.prompter.dialogs.back().find("detail")->stringOr(),
             std::string("Input & cleaning › Frames › Stop: 'x' is not a whole number"));
}

void WindowTests::undoRedo() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(!w.action(QStringLiteral("undo"))->isEnabled());
    w.document().set("analyses.contacts.cutoff", 7.5);
    QVERIFY(w.action(QStringLiteral("undo"))->isEnabled());
    QTest::keyClick(&w, Qt::Key_Z, Qt::ControlModifier);
    QCOMPARE(w.document().text("analyses.contacts.cutoff"), std::string("4.5"));
    QVERIFY(w.action(QStringLiteral("redo"))->isEnabled());
    w.action(QStringLiteral("redo"))->trigger();
    QCOMPARE(w.document().text("analyses.contacts.cutoff"), std::string("7.5"));
}

void WindowTests::projectLog() {
    Fixture f;
    app::MainWindow& w = *f.window;
    w.action(QStringLiteral("project_log"))->trigger();
    QCOMPARE(f.prompter.dialogs.back().find("text")->stringOr(), std::string("The project has no file yet"));
    QVERIFY(w.openPath(f.copy("default.mdatk.toml", "p.mdatk.toml")));
    w.action(QStringLiteral("project_log"))->trigger();
    QVERIFY(w.lastDialog());
    QCOMPARE(w.lastDialog()->windowTitle(), QStringLiteral("Project log"));
    // a project with a step log: MDATK's lines (time  step  [status]  title)
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    w.action(QStringLiteral("project_log"))->trigger();
    const QString text = w.lastDialog()->text();
    QVERIFY2(text.contains(QStringLiteral("  [")) || text == QStringLiteral("Nothing has been run in this project yet."),
             qPrintable(text));
    w.lastDialog()->close();
}

void WindowTests::citations() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "p.mdatk.toml")));
    const app::ProjectCitations expected = app::projectCitations(w.document().resolvedProject());
    w.action(QStringLiteral("citations"))->trigger();
    QVERIFY(w.lastDialog());
    QCOMPARE(w.lastDialog()->text().toStdString(), app::citationText(expected));
    w.lastDialog()->findChild<QPushButton*>(QStringLiteral("dialog.copy"))->click();
    QCOMPARE(QGuiApplication::clipboard()->text().toStdString(), app::citationText(expected));
    w.lastDialog()->findChild<QPushButton*>(QStringLiteral("dialog.copy_json"))->click();
    QVERIFY(QGuiApplication::clipboard()->text().contains(QStringLiteral("\"settings_keys\"")));
    w.lastDialog()->close();
    // Export citations (BibTeX)
    const fs::path bib = f.dir / "refs.bib";
    f.nextFile = {qs(bib)};
    w.action(QStringLiteral("export_citations"))->trigger();
    if (expected.keys.empty()) {
        QCOMPARE(f.prompter.dialogs.back().find("text")->stringOr(), std::string("Nothing to cite yet"));
    } else {
        QVERIFY(f.requests.back().start.endsWith(QStringLiteral("_citations.bib")));
        QCOMPARE(readFile(bib), core::citations::bibtex(expected.keys));
        QCOMPARE(f.lastLog(), "Wrote " + std::to_string(expected.keys.size()) + " references to " + shownText(bib));
    }
    // a project without any resource: nothing to cite
    w.action(QStringLiteral("new"))->trigger();
    const app::ProjectCitations none = app::projectCitations(w.document().resolvedProject());
    if (none.keys.empty()) {
        w.action(QStringLiteral("export_citations"))->trigger();
        QCOMPARE(f.prompter.dialogs.back().find("text")->stringOr(), std::string("Nothing to cite yet"));
    }
}

void WindowTests::aboutAndGraphics() {
    Fixture f;
    app::MainWindow& w = *f.window;
    w.action(QStringLiteral("about"))->trigger();
    QVERIFY(w.lastDialog());
    QCOMPARE(w.lastDialog()->windowTitle(), QStringLiteral("About MDBridge"));
    QVERIFY(w.lastDialog()->findChild<QLabel*>(QStringLiteral("dialog.heading"))->text().contains(QString::fromLatin1(mdbridge::kVersion)));
    QVERIFY(w.lastDialog()->findChild<QPushButton*>(QStringLiteral("dialog.bibtex")));
    w.action(QStringLiteral("graphics"))->trigger();
    QCOMPARE(w.lastDialog()->windowTitle(), QStringLiteral("Graphics"));
    QVERIFY(w.lastDialog()->text().contains(QStringLiteral("3D renderer chosen: software")));
    w.lastDialog()->close();
}

void WindowTests::workflowShells() {
    Fixture f;
    // D2: present, disabled, saying when they come (milestone 7)
    for (const char* name : {"planner", "wizard"}) {
        QAction* a = f.window->action(QString::fromLatin1(name));
        QVERIFY2(a && !a->isEnabled(), name);
        QVERIFY(a->toolTip().contains(QStringLiteral("milestone 7")));
    }
    // milestone 6 W3-B (D1): the PDB clean-up and the structure builder work (tests/gui/test_tools.cpp)
    for (const char* name : {"pdb_tools", "builder"}) {
        QAction* a = f.window->action(QString::fromLatin1(name));
        QVERIFY2(a && a->isEnabled(), name);
    }
}

void WindowTests::demoProject() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path folder = f.dir / "demo";
    f.nextFile = {qs(folder)};
    w.action(QStringLiteral("demo_project"))->trigger();
    QCOMPARE(f.requests.back().mode, app::forms::PathMode::Folder);
    QVERIFY(w.jobRunner()->running() || fs::exists(folder / "demo.mdatk.toml"));
    QTRY_VERIFY_WITH_TIMEOUT(!w.jobRunner()->running(), 300000);
    QTRY_COMPARE(w.document().path(), folder / "demo.mdatk.toml");
    const std::vector<std::string> log = f.log();
    QVERIFY(std::find(log.begin(), log.end(), "--- Creating demo ---") != log.end());
    QVERIFY(std::find(log.begin(), log.end(), "Demo system ready: " + shownText(folder / "demo.mdatk.toml")) !=
            log.end());
    QCOMPARE(f.lastLog(), "Opened project " + shownText(folder / "demo.mdatk.toml"));
}

// AP-12 (MDATK 0.1.0): the start with a file that does not load says so in the window
void WindowTests::startWithBrokenFile() {
    Fixture f;
    app::MainWindow& w = *f.window;
    w.logView()->setExpanded(false);
    const fs::path broken = f.copy("broken.mdatk.toml", "broken.mdatk.toml");
    QVERIFY(!w.startWith(broken));
    QCOMPARE(f.prompter.stderrLines.size(), std::size_t{1});
    QVERIFY(f.prompter.stderrLines[0].rfind("warning: could not load " + core::pathToUtf8(broken) + ": ", 0) == 0);
    const std::vector<std::string> log = f.log();
    QCOMPARE(log.size(), std::size_t{2});
    QVERIFY(log[0].rfind("ERROR: Could not open " + shownText(broken) + ": ", 0) == 0);
    QCOMPARE(log[1], std::string("These are the default settings instead (unsaved)."));
    QVERIFY(w.logView()->expanded());
    QVERIFY(f.prompter.dialogs.empty());
    QVERIFY(!w.document().hasPath());
    // a good file: no log line, not in the recent list (MDATK's start)
    QVERIFY(w.startWith(f.copy("default.mdatk.toml", "p.mdatk.toml")));
    QCOMPARE(f.log().size(), std::size_t{2});
}

// AP-13: closing with unsaved changes asks
void WindowTests::closingAsks() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("default.mdatk.toml", "p.mdatk.toml")));
    w.document().set("notes", std::string("unsaved"));
    f.answer("unsaved", "cancel");
    QVERIFY(!w.close());
    QVERIFY(w.isVisible());
    f.answer("unsaved", "discard");
    QVERIFY(w.close());
    QVERIFY(!w.isVisible());
    QVERIFY(readFile(f.dir / "p.mdatk.toml").find("unsaved") == std::string::npos);
}

// AP-14: closing during a job asks; Close cancels the job and waits for its clean-up
void WindowTests::closingDuringAJob() {
    Fixture f;
    app::MainWindow& w = *f.window;
    auto started = std::make_shared<std::atomic<bool>>(false);
    auto cleaned = std::make_shared<std::atomic<bool>>(false);
    QVERIFY(w.jobRunner()->start(
        QStringLiteral("Waiting"),
        [started, cleaned](const core::Reporter& reporter) {
            started->store(true);
            try {
                for (;;) {
                    reporter.progress(0.5, "waiting");
                    reporter.checkCancel();
                    std::this_thread::sleep_for(std::chrono::milliseconds(5));
                }
            } catch (const core::Cancelled&) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));  // the job's own clean-up
                cleaned->store(true);
                throw;
            }
        },
        [] {}));
    QVERIFY(w.action(QStringLiteral("cancel"))->isEnabled());
    // a job cancelled before its thread picked it up never runs (no clean-up to wait for): let it start
    QTRY_VERIFY_WITH_TIMEOUT(started->load(), 10000);
    int asked = 0;
    w.jobRunner()->setAsker([&asked](QWidget*, app::jobs::CloseQuestion, const QString&) {
        ++asked;
        return asked == 1 ? app::jobs::CloseAnswer::Cancel : app::jobs::CloseAnswer::Close;
    });
    QVERIFY(!w.close());
    QVERIFY(w.jobRunner()->running());
    const bool closed = w.close();
    QVERIFY(closed);
    QVERIFY(cleaned->load());
    QVERIFY(!w.jobRunner()->running());
}

// D6: a recovery copy of a modified project, offered at the next start
void WindowTests::recoveryCopies() {
    fs::path recovery;
    fs::path project;
    {
        Fixture f;
        app::MainWindow& w = *f.window;
        project = f.copy("default.mdatk.toml", "p.mdatk.toml");
        fs::copy_file(project, core::pathFromUtf8(m_config.path().toStdString()) / "p.mdatk.toml",
                      fs::copy_options::overwrite_existing);
        project = core::pathFromUtf8(m_config.path().toStdString()) / "p.mdatk.toml";
        QVERIFY(w.openPath(project));
        w.writeRecoveryNow();
        QVERIFY(w.recoveryFile().empty());  // unmodified: none
        w.document().set("analyses.contacts.cutoff", 9.25);
        w.writeRecoveryNow();
        recovery = w.recoveryFile();
        QVERIFY(fs::exists(recovery));
        QCOMPARE(recovery.parent_path(), core::recoveryDir());
        // the window "crashes": the copy stays
        w.setPrompter(nullptr);
        f.window->cancelLoad();
        f.window.reset();
    }
    QVERIFY(fs::exists(recovery));
    // a copy written moments ago may be another open window's: not offered
    {
        Fixture f;
        f.window->offerRecoveries();
        QVERIFY(f.prompter.questions.empty());
    }
    fs::last_write_time(recovery, fs::file_time_type::clock::now() - std::chrono::minutes(10));
    {
        Fixture f;
        f.answer("recover", "recover");
        f.window->offerRecoveries();
        QCOMPARE(f.lastQuestion(), std::string("recover:recover"));
        QCOMPARE(f.window->document().path(), project);
        QVERIFY(f.window->document().isModified());
        QCOMPARE(f.window->document().text("analyses.contacts.cutoff"), std::string("9.25"));
        QVERIFY(!fs::exists(recovery) || f.window->recoveryFile() == recovery);  // replaced by the window's own
        // saving drops the copy
        QVERIFY(f.window->save());
        QVERIFY(f.window->recoveryFile().empty());
        QVERIFY(core::findRecoveries().empty());
    }
}

void WindowTests::externalChange() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path p = f.copy("default.mdatk.toml", "p.mdatk.toml");
    QVERIFY(w.openPath(p));
    QVERIFY(!w.externalChangeShown());
    {
        std::ofstream out(p, std::ios::app | std::ios::binary);
        out << "\n# changed by another program\n";
    }
    w.checkDisk();
    QVERIFY(w.externalChangeShown());
    QCOMPARE(QString::fromStdString(f.lastLog()), shown(p) + QStringLiteral(" was changed by another program."));
    // Save asks before it replaces the other program's version
    w.document().set("notes", std::string("mine"));
    f.answer("external", "cancel");
    QVERIFY(!w.save());
    QCOMPARE(f.lastQuestion(), std::string("external:cancel"));
    // Reload (discarding the edit) reads the file again
    f.answer("reload", "reload");
    w.findChild<QPushButton*>(QStringLiteral("external.reload"))->click();
    QVERIFY(!w.document().isModified());
    QVERIFY(!w.document().changedOnDisk());
    QVERIFY(!w.externalChangeShown());
}

void WindowTests::dragAndDrop() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path p = f.copy("default.mdatk.toml", "dropped.mdatk.toml");
    QMimeData data;
    data.setUrls({QUrl::fromLocalFile(qs(p))});
    QDragEnterEvent enter(QPoint(10, 10), Qt::CopyAction, &data, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&w, &enter);
    QVERIFY(enter.isAccepted());
    QDropEvent drop(QPointF(10, 10), Qt::CopyAction, &data, Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&w, &drop);
    QTRY_COMPARE(w.document().path(), p);
}

// Every page fits the side bar (WN-11): with all sections open, no widget of a page reaches past
// its scroll area's viewport (rows wrap, hints wrap, tables fit their columns) and the tables
// need no horizontal scroll bar.
void WindowTests::pagesFitTheSideBar() {
    Fixture f;
    app::MainWindow& w = *f.window;
    // measured in the application's font, the bundled IBM Plex Serif (FONT-01; wider than DejaVu Sans)
    QCOMPARE(QFontInfo(w.font()).family(), QStringLiteral("IBM Plex Serif"));
    // D10: also run at 150 % and 200 % (ctest gui.window.scale150 / scale200: QT_SCALE_FACTOR)
    if (const QByteArray factor = qgetenv("QT_SCALE_FACTOR"); !factor.isEmpty()) {
        QCOMPARE(w.devicePixelRatioF(), factor.toDouble());
    }
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    QStringList overflows;
    for (const app::MainWindow::Page& page : w.pages()) {
        QVERIFY(w.showPage(page.id));
        for (QToolButton* header : page.widget->findChildren<QToolButton*>()) {
            if (header->objectName().startsWith(QStringLiteral("section:")) && header->isCheckable() &&
                !header->isChecked()) {
                header->setChecked(true);
            }
        }
        QTest::qWait(50);  // the layouts settle (tables fit their columns after a reload)
        auto* scroll = page.widget->findChild<QScrollArea*>(QStringLiteral("page:%1:scroll").arg(page.id));
        if (!scroll) {
            QVERIFY2(page.placeholder, qPrintable(page.id));
            continue;
        }
        QWidget* viewport = scroll->viewport();
        QWidget* body = scroll->widget();
        QVERIFY(body);
        const int room = viewport->width();
        QVERIFY2(room > 400, qPrintable(QStringLiteral("%1: viewport %2 px").arg(page.id).arg(room)));
        if (body->width() > room) {
            overflows << QStringLiteral("%1: the page is %2 px wide in a %3 px viewport (minimum %4 px)")
                             .arg(page.id)
                             .arg(body->width())
                             .arg(room)
                             .arg(body->minimumSizeHint().width());
        }
        int checked = 0;
        for (QWidget* child : body->findChildren<QWidget*>()) {
            if (!child->isVisibleTo(body)) {
                continue;
            }
            bool insideView = false;  // a table's cells and scroll bars are the table's business
            for (QWidget* up = child->parentWidget(); up && up != body; up = up->parentWidget()) {
                insideView = insideView || qobject_cast<QAbstractScrollArea*>(up) != nullptr;
            }
            if (insideView) {
                continue;
            }
            ++checked;
            // no text cut off: buttons and labels at least their minimum size (also at 150 % and 200 %)
            const auto* label = qobject_cast<QLabel*>(child);
            const bool own = !child->inherits("QLineEditIconButton");  // Qt sizes a line edit's icon itself
            if (own && (qobject_cast<QAbstractButton*>(child) || (label && !label->wordWrap() && !label->text().isEmpty())) &&
                (child->width() < child->minimumSizeHint().width() || child->height() < child->minimumSizeHint().height())) {
                overflows << QStringLiteral("%1: %2 '%3' is %4 x %5 px, less than its %6 x %7 px")
                                 .arg(page.id, QString::fromLatin1(child->metaObject()->className()), child->objectName())
                                 .arg(child->width())
                                 .arg(child->height())
                                 .arg(child->minimumSizeHint().width())
                                 .arg(child->minimumSizeHint().height());
            }
            const int right = child->mapTo(viewport, QPoint(child->width(), 0)).x();
            if (right > room) {
                overflows << QStringLiteral("%1: %2 '%3' ends at %4 px (minimum width %5 px)")
                                 .arg(page.id, QString::fromLatin1(child->metaObject()->className()), child->objectName())
                                 .arg(right)
                                 .arg(child->minimumSizeHint().width());
            }
            if (auto* table = qobject_cast<QTableView*>(child); table && table->horizontalScrollBar()->isVisible()) {
                overflows << QStringLiteral("%1: the table '%2' scrolls sideways").arg(page.id, table->objectName());
            }
            if (auto* table = qobject_cast<QTableView*>(child)) {  // no column label clipped
                const QHeaderView* header = table->horizontalHeader();
                for (int c = 0; c < header->count(); ++c) {
                    if (!header->isSectionHidden(c) && header->sectionSize(c) < header->sectionSizeHint(c)) {
                        overflows << QStringLiteral("%1: the table '%2' clips its column label %3 (%4 of %5 px)")
                                         .arg(page.id, table->objectName())
                                         .arg(c)
                                         .arg(header->sectionSize(c))
                                         .arg(header->sectionSizeHint(c));
                    }
                }
            }
        }
        QVERIFY2(checked > 10, qPrintable(page.id));
    }
    // the window's own bars at MDATK's 1320 px: no button or label cut off
    for (const char* part : {"header", "pagebar", "actionbar"}) {
        for (QWidget* child : w.findChild<QWidget*>(QString::fromLatin1(part))->findChildren<QWidget*>()) {
            const auto* label = qobject_cast<QLabel*>(child);
            if (!child->isVisibleTo(&w) || !(qobject_cast<QAbstractButton*>(child) || (label && !label->text().isEmpty()))) continue;
            if (child->width() < child->minimumSizeHint().width() || child->height() < child->minimumSizeHint().height()) {
                overflows << QStringLiteral("%1: %2 '%3' is %4 x %5 px, less than its %6 x %7 px")
                                 .arg(QString::fromLatin1(part), QString::fromLatin1(child->metaObject()->className()),
                                      child->objectName())
                                 .arg(child->width())
                                 .arg(child->height())
                                 .arg(child->minimumSizeHint().width())
                                 .arg(child->minimumSizeHint().height());
            }
        }
    }
    QVERIFY2(overflows.isEmpty(), qPrintable(overflows.join(QLatin1Char('\n'))));
}

// WN-12: paths in interface text with the platform's separators. The Windows rules, checked on every
// platform; on Linux and macOS the window's helpers return their input.
void WindowTests::pathsShownNatively() {
    using app::gui::PathStyle;
    const auto win = [](const QString& text, const QStringList& roots = {}) {
        return app::gui::displayPathsInText(text, PathStyle::Windows, roots);
    };
    QCOMPARE(app::gui::displayPathText(QStringLiteral("C:/Users/me/out\\video.mp4"), PathStyle::Windows),
             QStringLiteral("C:\\Users\\me\\out\\video.mp4"));
    QCOMPARE(app::gui::displayPathText(QStringLiteral("C:/Users/me/out\\video.mp4"), PathStyle::Posix),
             QStringLiteral("C:/Users/me/out\\video.mp4"));
    // the pipeline's mixed joins, a path ending at ':' or white space, text around it unchanged
    QCOMPARE(win(QStringLiteral("MDBridge 0.5 - output directory C:/Users/me/out")),
             QStringLiteral("MDBridge 0.5 - output directory C:\\Users\\me\\out"));
    QCOMPARE(win(QStringLiteral("Video: C:/Users/me/out\\aptamer.mp4 (120/120 frames)")),
             QStringLiteral("Video: C:\\Users\\me\\out\\aptamer.mp4 (120/120 frames)"));
    QCOMPARE(win(QStringLiteral("ERROR: Could not open C:/p/a.mdatk.toml: no such file")),
             QStringLiteral("ERROR: Could not open C:\\p\\a.mdatk.toml: no such file"));
    QCOMPARE(win(QStringLiteral("Density map not found: d:/maps/x.mrc (untick 'show' to render without it)")),
             QStringLiteral("Density map not found: d:\\maps\\x.mrc (untick 'show' to render without it)"));
    QCOMPARE(win(QStringLiteral("skipping data source 'C:/x/e.dat': bad")), QStringLiteral("skipping data source 'C:\\x\\e.dat': bad"));
    // UNC and device paths; no path: unchanged
    QCOMPARE(win(QStringLiteral("Opened //server/share/p/x.nc")), QStringLiteral("Opened \\\\server\\share\\p\\x.nc"));
    QCOMPARE(win(QStringLiteral("at \\\\?\\C:/long/x.nc")), QStringLiteral("at \\\\?\\C:\\long\\x.nc"));
    for (const char* plain : {"5/10 frames, kcal/mol", "see https://example.org/a/b", "group A/B: 3 atoms", "inputs/a.parm7",
                              "AB:/x"}) {
        QCOMPARE(win(QString::fromLatin1(plain)), QString::fromLatin1(plain));
    }
    // a root folder with spaces is converted as a whole (the home folder, the project's folders)
    QCOMPARE(win(QStringLiteral("Wrote 3 references to C:/Users/John Doe/p/x.bib"), {QStringLiteral("C:\\Users\\John Doe")}),
             QStringLiteral("Wrote 3 references to C:\\Users\\John Doe\\p\\x.bib"));
    QCOMPARE(win(QStringLiteral("Saved c:/users/john doe/p.mdatk.toml"), {QStringLiteral("C:/Users/John Doe/")}),
             QStringLiteral("Saved c:\\users\\john doe\\p.mdatk.toml"));
    QCOMPARE(win(QStringLiteral("C:/database/a b/c"), {QStringLiteral("C:/data")}), QStringLiteral("C:\\database\\a b/c"));
    // the log view shows its lines through the window's transform
    app::jobs::LogView view;
    view.setLineTransform([&](const QString& line) { return win(line); });
    view.feed(QStringLiteral("Opened C:/a\\b.mp4"));
    QCOMPARE(view.text(), QStringLiteral("Opened C:\\a\\b.mp4\n"));
    if constexpr (app::gui::kNativePathStyle == PathStyle::Posix) {
        QCOMPARE(app::gui::displayPath(fs::path("/home/me/a b/p.mdatk.toml")), QStringLiteral("/home/me/a b/p.mdatk.toml"));
        QCOMPARE(app::gui::displayPathsIn(QStringLiteral("C:/x")), QStringLiteral("C:/x"));
    } else {
        QCOMPARE(app::gui::displayPath(fs::path("C:/me/out") / "v.mp4"), QStringLiteral("C:\\me\\out\\v.mp4"));
    }
    // the window: title, recent list, log, status bar
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path p = f.copy("default.mdatk.toml", "p.mdatk.toml");
    QVERIFY(w.openPath(p));
    QVERIFY(w.windowTitle().startsWith(shown(p)));
    QCOMPARE(QString::fromStdString(f.lastLog()), QStringLiteral("Opened project ") + shown(p));
    QCOMPARE(w.statusBar()->currentMessage(), QStringLiteral("Opened ") + shown(p));
    if constexpr (app::gui::kNativePathStyle == PathStyle::Windows) {
        QVERIFY(!w.windowTitle().contains(QLatin1Char('/')));
        QVERIFY(!QString::fromStdString(f.lastLog()).contains(QLatin1Char('/')));
    }
}

namespace {

// The accessible name a screen reader announces for `w` (its own, its text, a buddy label's).
QString accessibleNameOf(QWidget* w) {
    QAccessibleInterface* iface = QAccessible::queryAccessibleInterface(w);
    return iface ? iface->text(QAccessible::Name).trimmed() : QString();
}

// Opens every section of every page (all fields shown).
void openAllSections(app::MainWindow& w) {
    for (const app::MainWindow::Page& page : w.pages()) {
        w.showPage(page.id);
        for (QToolButton* header : page.widget->findChildren<QToolButton*>()) {
            if (header->objectName().startsWith(QStringLiteral("section:")) && header->isCheckable() && !header->isChecked()) {
                header->setChecked(true);
            }
        }
    }
    QApplication::processEvents();
}

// Where a widget is: its object name and those of its named parents.
QString whereIs(const QWidget* w) {
    QStringList names;
    for (const QWidget* up = w; up; up = up->parentWidget()) {
        if (!up->objectName().isEmpty()) names.push_front(up->objectName());
        if (names.size() >= 3) break;
    }
    return QStringLiteral("%1 '%2'").arg(QString::fromLatin1(w->metaObject()->className()), names.join(QStringLiteral(" > ")));
}

// The widget the focus goes to for `w` (a spin box's line edit hands it to the spin box).
QWidget* focusOwner(QWidget* w) {
    while (w->focusProxy()) w = w->focusProxy();
    return w;
}

}  // namespace

// D10: every widget that takes the keyboard focus, on every page with every section open and in the
// window's chrome, has an accessible name (what a screen reader announces: its own name, its text,
// or a buddy label's).
void WindowTests::everyFocusableWidgetIsNamed() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    QStringList unnamed;
    int checked = 0;
    const auto check = [&](QWidget* root) {
        for (QWidget* child : root->findChildren<QWidget*>()) {
            // a widget with a focus proxy hands the focus on (a spin box's line edit to the spin box)
            if (!child->isVisibleTo(&w) || !(child->focusPolicy() & Qt::TabFocus) || child->focusProxy()) continue;
            ++checked;
            if (accessibleNameOf(child).isEmpty()) unnamed << whereIs(child);
        }
    };
    openAllSections(w);
    for (const app::MainWindow::Page& page : w.pages()) {
        w.showPage(page.id);
        QApplication::processEvents();
        check(page.widget);
    }
    w.logView()->setExpanded(true);
    QApplication::processEvents();
    for (const char* part : {"header", "pagebar", "preview", "log", "actionbar"}) {
        QWidget* widget = w.findChild<QWidget*>(QString::fromLatin1(part));
        QVERIFY2(widget, part);
        if (widget->focusPolicy() & Qt::TabFocus) {
            ++checked;
            if (accessibleNameOf(widget).isEmpty()) unnamed << whereIs(widget);
        }
        check(widget);
    }
    QVERIFY2(checked > 300, qPrintable(QString::number(checked)));
    unnamed.removeDuplicates();
    QVERIFY2(unnamed.isEmpty(), qPrintable(unnamed.join(QLatin1Char('\n'))));
}

// D10: every control (button, field, list, slider, the canvas) can take the keyboard focus with
// Tab: none is reachable by the mouse only. The header's tool buttons too (a QToolBar's have no
// focus by default).
void WindowTests::everyControlTakesTheKeyboardFocus() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    openAllSections(w);
    w.logView()->setExpanded(true);
    QStringList unreachable;
    int controls = 0;
    QStringList seen;
    for (const app::MainWindow::Page& page : w.pages()) {
        w.showPage(page.id);
        QApplication::processEvents();
        for (QWidget* child : w.findChildren<QWidget*>()) {
            if (!child->isVisibleTo(&w)) continue;
            const bool control = qobject_cast<QAbstractButton*>(child) || qobject_cast<QComboBox*>(child) ||
                                 qobject_cast<QLineEdit*>(child) || qobject_cast<QAbstractSpinBox*>(child) ||
                                 qobject_cast<QAbstractSlider*>(child) || qobject_cast<QAbstractItemView*>(child) ||
                                 qobject_cast<QPlainTextEdit*>(child) || qobject_cast<QTextEdit*>(child) ||
                                 qobject_cast<app::PreviewWidget*>(child);
            // not controls of their own: scroll bars, a table's header, a line edit's icon (its action)
            if (!control || qobject_cast<QScrollBar*>(child) || qobject_cast<QHeaderView*>(child) ||
                child->inherits("QLineEditIconButton")) {
                continue;
            }
            if (child->objectName().startsWith(QStringLiteral("qt_"))) continue;  // Qt's own parts (a line edit's clear button)
            const QString where = whereIs(child);
            if (seen.contains(where)) continue;
            seen << where;
            ++controls;
            QWidget* owner = focusOwner(child);
            const bool takes = (owner->focusPolicy() & Qt::TabFocus) || (child->focusPolicy() & Qt::TabFocus) ||
                               (child->focusProxy() && (child->focusProxy()->focusPolicy() & Qt::TabFocus));
            if (!takes) unreachable << where;
        }
    }
    QVERIFY2(controls > 300, qPrintable(QString::number(controls)));
    QVERIFY2(unreachable.isEmpty(), qPrintable(unreachable.join(QLatin1Char('\n'))));
    // Tab goes round the window: header, page buttons, the page, the canvas, the log, the action bar
    w.showPage(QStringLiteral("groups"));
    QApplication::processEvents();
    QWidget* start = w.findChild<QWidget*>(QStringLiteral("header.open"));
    QVERIFY(start);
    w.activateWindow();
    start->setFocus(Qt::TabFocusReason);
    QTRY_COMPARE(QApplication::focusWidget(), start);
    QStringList order;
    for (int step = 0; step < 3000; ++step) {
        QWidget* now = QApplication::focusWidget();
        QVERIFY(now);
        order << now->objectName();
        QTest::keyClick(now, Qt::Key_Tab);
        if (QApplication::focusWidget() == start) break;
    }
    QCOMPARE(QApplication::focusWidget(), start);
    QStringList expected{QStringLiteral("header.open"), QStringLiteral("header.workflow"), QStringLiteral("header.menu"),
                         QStringLiteral("page:input"), QStringLiteral("pagebar.refresh"), QStringLiteral("groups"),
                         QStringLiteral("preview"), QStringLiteral("log.text"), QStringLiteral("actionbar.render"),
                         QStringLiteral("actionbar.open_video")};
    qsizetype at = -1;
    for (const QString& name : expected) {
        const qsizetype next = order.indexOf(name);
        QVERIFY2(next > at, qPrintable(name + QStringLiteral(" in ") + order.join(QStringLiteral(", "))));
        at = next;
    }
    // the page buttons are one group: Tab enters it, the arrow keys move in it, Space shows the page
    QWidget* input = w.findChild<QWidget*>(QStringLiteral("page:input"));
    input->setFocus(Qt::TabFocusReason);
    QTRY_COMPARE(QApplication::focusWidget(), input);
    QTest::keyClick(input, Qt::Key_Right);
    QCOMPARE(QApplication::focusWidget()->objectName(), QStringLiteral("page:groups"));
    QTest::keyClick(QApplication::focusWidget(), Qt::Key_Space);
    QCOMPARE(w.currentPage(), QStringLiteral("groups"));
}

// D10: on every page, Tab from the page's first field visits every field editor shown (every section
// open), section by section in the page's order, without being caught by a field (a multi-line text
// passes the focus on; a table's cells are left at its end).
void WindowTests::tabOrderVisitsEveryField() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    openAllSections(w);
    QStringList problems;
    for (const app::MainWindow::Page& page : w.pages()) {
        if (page.placeholder) continue;
        w.showPage(page.id);
        QApplication::processEvents();
        std::vector<QWidget*> editors;
        for (app::forms::FieldWidget* field : page.widget->findChildren<app::forms::FieldWidget*>()) {
            QWidget* e = field->editor();
            if (e && e->isVisibleTo(&w) && e->isEnabled()) editors.push_back(e);
        }
        QVERIFY2(!editors.empty(), qPrintable(page.id));
        // the sections in the page's order
        QList<QWidget*> sections;
        for (QWidget* box : page.widget->findChildren<QWidget*>()) {
            if (box->objectName().startsWith(QStringLiteral("section:")) && !qobject_cast<QToolButton*>(box)) sections << box;
        }
        const auto sectionOf = [&](QWidget* e) {
            for (qsizetype k = 0; k < sections.size(); ++k) {
                if (sections[k]->isAncestorOf(e)) return static_cast<int>(k);
            }
            return -1;
        };
        // the field editor a focused widget belongs to (a spin box's line edit, a table's cell editor)
        const auto editorOf = [&](QWidget* c) -> QWidget* {
            for (QWidget* up = c; up && up != page.widget; up = up->parentWidget()) {
                if (std::find(editors.begin(), editors.end(), up) != editors.end()) return up;
            }
            return nullptr;
        };
        // start at the first field in the focus chain of the page
        QWidget* first = nullptr;
        for (QWidget* c = page.widget->nextInFocusChain(); c && c != page.widget; c = c->nextInFocusChain()) {
            if (page.widget->isAncestorOf(c) && editorOf(c) && (c->focusPolicy() & Qt::TabFocus) && c->isVisibleTo(&w)) {
                first = c;
                break;
            }
        }
        QVERIFY2(first, qPrintable(page.id));
        w.activateWindow();
        first->setFocus(Qt::TabFocusReason);
        QTRY_COMPARE(QApplication::focusWidget(), first);
        std::vector<QWidget*> visited;
        int lastSection = -1;
        for (int step = 0; step < 4000; ++step) {
            QWidget* now = QApplication::focusWidget();
            if (!now || !page.widget->isAncestorOf(now)) break;
            QWidget* owner = editorOf(now);
            if (owner && std::find(visited.begin(), visited.end(), owner) == visited.end()) {
                visited.push_back(owner);
                const int section = sectionOf(owner);
                if (section >= 0 && section < lastSection) {
                    problems << QStringLiteral("%1: Tab goes back from section %2 to %3 at %4")
                                    .arg(page.id)
                                    .arg(lastSection)
                                    .arg(section)
                                    .arg(owner->objectName());
                }
                lastSection = std::max(lastSection, section);
            }
            QTest::keyClick(now, Qt::Key_Tab);
        }
        if (QWidget* now = QApplication::focusWidget(); now && page.widget->isAncestorOf(now)) {
            problems << QStringLiteral("%1: Tab is caught at %2").arg(page.id, whereIs(now));
        }
        for (QWidget* e : editors) {
            if (std::find(visited.begin(), visited.end(), e) == visited.end()) {
                problems << QStringLiteral("%1: Tab does not reach %2").arg(page.id, e->objectName());
            }
        }
    }
    QVERIFY2(problems.isEmpty(), qPrintable(problems.join(QLatin1Char('\n'))));
}

// D10: Esc closes the window's dialogs (About, Citations, Graphics, Project log).
void WindowTests::escapeClosesDialogs() {
    Fixture f;
    app::MainWindow& w = *f.window;
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));  // the project log needs a file
    for (const char* name : {"about", "citations", "graphics", "project_log"}) {
        w.action(QString::fromLatin1(name))->trigger();
        QDialog* dialog = nullptr;
        QTRY_VERIFY((dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget())) != nullptr ||
                    [&] {
                        for (QDialog* d : w.findChildren<QDialog*>()) {
                            if (d->isVisible()) dialog = d;
                        }
                        return dialog != nullptr;
                    }());
        QVERIFY2(dialog, name);
        QTest::keyClick(dialog, Qt::Key_Escape);
        QTRY_VERIFY2(!dialog->isVisible(), name);
    }
}

// D16: the recent list in MDATK's shared settings file on Linux ($XDG_CONFIG_HOME/mdatk), MDBridge's
// own settings and recovery copies in $XDG_CONFIG_HOME/mdbridge; on Windows both in
// %APPDATA%\MDBridge (MDATK's Windows build keeps its own file elsewhere, not shared).
void WindowTests::settingsLocations() {
    const fs::path config = core::pathFromUtf8(m_config.path().toStdString());
#ifdef _WIN32
    const fs::path recentFile = config / "MDBridge" / "settings.toml";
    const fs::path ownDir = config / "MDBridge";
#else
    const fs::path recentFile = config / "mdatk" / "settings.toml";
    const fs::path ownDir = config / "mdbridge";
#endif
    QCOMPARE(core::userSettingsPath(), recentFile);
    QCOMPARE(core::mdbridgeConfigDir(), ownDir);
    QCOMPARE(core::recoveryDir(), ownDir / "recovery");
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path p = f.copy("default.mdatk.toml", "where.mdatk.toml");
    QVERIFY(w.openPath(p));
    QVERIFY(fs::exists(recentFile));
    QVERIFY(readFile(recentFile).find("where.mdatk.toml") != std::string::npos);
    w.document().set("notes", std::string("unsaved"));
    w.writeRecoveryNow();
    QVERIFY(!w.recoveryFile().empty());
    QCOMPARE(w.recoveryFile().parent_path(), ownDir / "recovery");
#ifndef _WIN32
    QVERIFY(!fs::exists(config / "MDBridge"));  // nothing in the Windows folder name
    QVERIFY(!fs::exists(config / "mdatk" / "recovery"));  // MDATK's folder holds the shared list only
#endif
    QVERIFY(w.save());
    QVERIFY(w.recoveryFile().empty());
}

void WindowTests::screenshots() {
    Fixture f;
    app::MainWindow& w = *f.window;
    const fs::path dir = guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR));
    QVERIFY(w.openPath(f.copy("full.mdatk.toml", "full.mdatk.toml")));
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading() && w.preview()->settled(), 120000);
    // at 150 % and 200 % (gui.window.scale150 / scale200) the names say so: window-input-150.png
    const QByteArray factor = qgetenv("QT_SCALE_FACTOR");
    const std::string suffix =
        factor.isEmpty() || factor == "1" ? std::string() : "-" + std::to_string(static_cast<int>(factor.toDouble() * 100.0 + 0.5));
    for (const app::MainWindow::Page& page : w.pages()) {
        w.showPage(page.id);
        QApplication::processEvents();
        const fs::path shot = guitest::saveScreenshot(&w, dir, "window-" + page.id.toStdString() + suffix);
        QVERIFY(fs::file_size(shot) > 10000);
        qInfo("screenshot %s", core::pathToUtf8(shot).c_str());
    }
}

int main(int argc, char** argv) {
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    QApplication::setStyle(QStringLiteral("Fusion"));
    guitest::useBundledFont();
    WindowTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_window.moc"
