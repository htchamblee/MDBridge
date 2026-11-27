// mdbridge_gui_tools_tests: the PDB clean-up and structure builder windows (src/app/gui/tools,
// milestone 6 W3-B, owner decision D1): opened from the analysis window's Workflow menu, hosted
// through tools::Host (also without the analysis window), the counterparts of MDATK's
// test_pdb_tools_window_checks_and_cleans and test_builder_non_wc_helix_builds_in_the_background,
// the files they write equal to what `mdbridge pdb clean` / `build helix` / `build mutate` / `build
// fit` write, Cancel of a fit logged as cancelled, closing during a build, the keyboard and
// accessibility rules of the analysis window (every control named and reachable with Tab, in the
// form's order), the forms fitting their column with nothing cut off (also at 150 % and 200 %:
// ctest gui.tools.scale150 / scale200) and screenshots (artifacts). The action scripts with MDATK's
// files are tests/data/gui/tools (mdbridge_gui_tests). Qt Test, offscreen.

#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <QAbstractButton>
#include <QAbstractItemView>
#include <QAbstractSpinBox>
#include <QAccessible>
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QSpinBox>
#include <QSplitter>
#include <QTabBar>
#include <QTabWidget>
#include <QTemporaryDir>
#include <QTest>
#include <QToolButton>

#include "commands.hpp"
#include "gui/dialogs/prompts.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "gui/main_window.hpp"
#include "gui/tools/builder_window.hpp"
#include "gui/tools/host.hpp"
#include "gui/tools/pdb_tools_window.hpp"
#include "gui/tools/tool_actions.hpp"
#include "gui/tools/tool_windows.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/io/pdb.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;
using app::tools::BuilderWindow;
using app::tools::PdbToolsWindow;
using app::tools::ToolWindow;

namespace {

const fs::path kData = core::pathFromUtf8(MDBRIDGE_TEST_DATA);
// MDATK's test_pdbtools make_pdb fixture of test_pdb_tools_window_checks_and_cleans (written by
// tests/parity/make_gui_cases.py): D G1 C2 (chain A), L 0G3 0U4 (chain B)
const fs::path kFourResidues = kData / "gui/tools/pdb_tools_check_and_clean/in/in.pdb";

fs::path shotDir() { return guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR)); }

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

QString qs(const fs::path& path) { return QString::fromStdString(core::pathToUtf8(path)); }

std::vector<std::string> residueNames(const std::string& text) {
    std::vector<std::string> out;
    const io::PdbStructure s = io::readPdbText(text);
    for (const io::PdbResidue* r : s.residues()) out.push_back(r->name);
    return out;
}

// `mdbridge <args>` in-process; its exit status
int runCommand(const std::vector<std::string>& args) {
    std::vector<const char*> argv{"mdbridge"};
    for (const std::string& a : args) argv.push_back(a.c_str());
    std::ostringstream out, err;
    app::Console console{out, err, false};
    const int status = app::runCoreCommandLine(static_cast<int>(argv.size()), argv.data(), console);
    if (status != 0) qWarning("mdbridge %s: %s", args.front().c_str(), err.str().c_str());
    return status;
}

// The step log entries of `project`, as tables.
core::Array stepLog(const fs::path& project) { return core::readLog(project); }

std::string entryText(const core::Value& entry, const char* key) {
    const core::Value* v = entry.find(key);
    return v ? v->stringOr() : std::string();
}

QString accessibleNameOf(QWidget* w) {
    QAccessibleInterface* iface = QAccessible::queryAccessibleInterface(w);
    return iface ? iface->text(QAccessible::Name).trimmed() : QString();
}

QString whereIs(const QWidget* w) {
    return QStringLiteral("%1 '%2'").arg(QString::fromLatin1(w->metaObject()->className()), w->objectName());
}

bool isControl(QWidget* w) {
    return qobject_cast<QAbstractButton*>(w) || qobject_cast<QComboBox*>(w) || qobject_cast<QLineEdit*>(w) ||
           qobject_cast<QAbstractSpinBox*>(w) || qobject_cast<QPlainTextEdit*>(w) || qobject_cast<QTabBar*>(w);
}

// A host that is not the analysis window (milestone 7's RNAdo will be one).
class FakeHost : public app::tools::Host {
public:
    QWidget* window() override { return nullptr; }
    core::ProjectDocument& document() override { return m_document; }
    fs::path projectPath() const override { return project; }
    void appendLog(const QString& text) override { log << text; }
    void error(const QString& text, const QString& detail) override { log << QStringLiteral("ERROR: ") + text + detail; }
    void info(const QString& text, const QString& detail) override { log << text + detail; }
    bool saveProject(bool) override { return false; }
    app::Prompter& prompter() override { return m_prompter; }

    fs::path project;
    QStringList log;

private:
    core::ProjectDocument m_document;
    guitest::RecordingPrompter m_prompter;
};

// The analysis window with a saved default project in a folder of its own, the current folder.
struct Fixture {
    QTemporaryDir tmp;
    fs::path dir;
    fs::path project;
    fs::path cwd;
    guitest::RecordingPrompter prompter;
    std::optional<app::MainWindow> window;

    Fixture() {
        dir = core::pathFromUtf8(tmp.path().toStdString());
        cwd = fs::current_path();
        fs::current_path(dir);
        project = dir / "p.mdatk.toml";
        fs::copy_file(kData / "projects/default.mdatk.toml", project);
        fs::copy_file(kFourResidues, dir / "in.pdb");
        window.emplace(render::Renderer::Software);
        window->setPrompter(&prompter);
        window->resize(1320, 920);
        window->show();
        if (!QTest::qWaitForWindowExposed(&*window)) qWarning("the window was not exposed");
        if (!window->openPath(project, false)) qWarning("the project did not open");
    }
    ~Fixture() {
        window->cancelLoad();
        window->setPrompter(nullptr);
        window.reset();
        fs::current_path(cwd);
    }
    app::tools::ToolWindows& tools() { return window->toolWindows(); }
    PdbToolsWindow* openPdbTools() {
        window->action(QStringLiteral("pdb_tools"))->trigger();
        return tools().pdbTools();
    }
    BuilderWindow* openBuilder() {
        window->action(QStringLiteral("builder"))->trigger();
        return tools().builder();
    }
    std::vector<std::string> log() {
        window->logView()->flush();
        std::vector<std::string> out;
        for (const QString& line : window->logView()->text().split(QLatin1Char('\n'), Qt::SkipEmptyParts)) {
            out.push_back(line.toStdString());
        }
        return out;
    }
};

template <class W>
W* control(ToolWindow* tool, const char* key) {
    auto* w = qobject_cast<W*>(tool->control(QString::fromLatin1(key)));
    if (!w) qWarning("no control %s", key);
    return w;
}

void setText(ToolWindow* tool, const char* key, const QString& text) { control<QLineEdit>(tool, key)->setText(text); }

void choose(ToolWindow* tool, const char* key, const char* id) {
    QComboBox* box = control<QComboBox>(tool, key);
    box->setCurrentIndex(box->findData(QString::fromLatin1(id)));
}

void click(ToolWindow* tool, const char* key) { QTest::mouseClick(control<QPushButton>(tool, key), Qt::LeftButton); }

// Every page of a tool window: the builder's three, the clean-up's one.
std::vector<int> pagesOf(ToolWindow* tool) {
    if (auto* builder = qobject_cast<BuilderWindow*>(tool)) {
        (void)builder;
        return {0, 1, 2};
    }
    return {0};
}

void showPage(ToolWindow* tool, int k) {
    if (auto* builder = qobject_cast<BuilderWindow*>(tool)) {
        builder->setPage(QString::fromLatin1(std::array{"helix", "mutate", "fit"}[static_cast<std::size_t>(k)]));
    }
    QApplication::processEvents();
}

QScrollArea* formScroll(ToolWindow* tool, int k) {
    const char* names[] = {"helix", "mutate", "fit"};
    const QString key = qobject_cast<BuilderWindow*>(tool) ? QString::fromLatin1(names[k]) : QStringLiteral("form");
    return tool->findChild<QScrollArea*>(tool->toolId() + QLatin1Char('.') + key + QStringLiteral(":scroll"));
}

}  // namespace

class ToolTests : public QObject {
    Q_OBJECT

private slots:
    void initTestCase();
    void menuOpensTheTools();
    void pdbToolsWindowChecksAndCleans();
    void pdbCleanWritesWhatTheCommandWrites();
    void builderHelixWritesWhatTheCommandWrites();
    void builderNonWcHelixBuildsInTheBackground();
    void mutateAndFitWriteWhatTheCommandsWrite();
    void fitCancelIsLoggedAsCancelled();
    void closingDuringABuildAsks();
    void hostedWithoutTheAnalysisWindow();
    void everyControlIsNamedAndReachable();
    void tabOrderFollowsTheForm();
    void formsFitTheirColumn();
    void screenshots();

private:
    QTemporaryDir m_config;
};

void ToolTests::initTestCase() {
    QVERIFY(m_config.isValid());
    qputenv("XDG_CONFIG_HOME", m_config.path().toUtf8());
    qputenv("APPDATA", m_config.path().toUtf8());
}

// D1: the menu entries open the windows (one of each, shown again), the planner and the wizard stay
// shells until milestone 7
void ToolTests::menuOpensTheTools() {
    Fixture f;
    for (const char* name : {"planner", "wizard"}) {
        QAction* a = f.window->action(QString::fromLatin1(name));
        QVERIFY2(a && !a->isEnabled(), name);
        QVERIFY(a->toolTip().contains(QStringLiteral("milestone 7")));
    }
    for (const char* name : {"pdb_tools", "builder"}) {
        QAction* a = f.window->action(QString::fromLatin1(name));
        QVERIFY2(a && a->isEnabled(), name);
        QVERIFY(!a->toolTip().contains(QStringLiteral("milestone")));
    }
    PdbToolsWindow* pdb = f.openPdbTools();
    QVERIFY(pdb && pdb->isVisible());
    QCOMPARE(pdb->windowTitle(), QStringLiteral("PDB clean-up"));
    QCOMPARE(pdb->parentWidget(), &*f.window);  // above the analysis window (MDATK's transient_for)
    QCOMPARE(f.openPdbTools(), pdb);            // the open window again
    BuilderWindow* builder = f.openBuilder();
    QVERIFY(builder && builder->isVisible());
    QCOMPARE(builder->windowTitle(), QStringLiteral("Structure builder"));
    QCOMPARE(builder->page(), QStringLiteral("helix"));
    pdb->close();
    QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    QVERIFY(!f.tools().pdbTools());  // MDATK: destroy sets window.pdb_tools to None
    QVERIFY(f.openPdbTools() != nullptr);
}

// MDATK tests/test_gui.py::test_pdb_tools_window_checks_and_cleans
void ToolTests::pdbToolsWindowChecksAndCleans() {
    Fixture f;
    PdbToolsWindow* dlg = f.openPdbTools();
    setText(dlg, "input", qs(f.dir / "in.pdb"));
    QVERIFY(control<QLineEdit>(dlg, "output")->text().endsWith(QStringLiteral("in_clean.pdb")));  // suggested
    dlg->runCheck();
    const QString report = dlg->reportText();
    QVERIFY2(report.contains(QStringLiteral("Chain B")) && report.contains(QStringLiteral("chirality 2 L")), qPrintable(report));
    choose(dlg, "scheme", "amber");
    control<QCheckBox>(dlg, "renumber")->setChecked(true);
    dlg->runClean();
    const fs::path out = core::pathFromUtf8(control<QLineEdit>(dlg, "output")->text().toStdString());
    const std::vector<std::string> names = residueNames(readFile(out));
    QCOMPARE(names, (std::vector<std::string>{"G5", "C3", "LG5", "LU3"}));
    QVERIFY(readFile(out).find("REMARK 950") == std::string::npos);
    QVERIFY(!control<QComboBox>(dlg, "ann_engine")->isEnabled());
    control<QCheckBox>(dlg, "annotate")->setChecked(true);
    choose(dlg, "ann_engine", "builtin");
    control<QCheckBox>(dlg, "ann_stacking")->setChecked(false);
    QVERIFY(control<QComboBox>(dlg, "ann_engine")->isEnabled());
    dlg->runClean();
    const std::string text = readFile(out);
    QVERIFY(text.find("REMARK 950 ENGINE builtin") != std::string::npos && text.find("REMARK 950 DBN") != std::string::npos);
    QCOMPARE(residueNames(text), names);
    // the host's log and the project's step log (MDATK's logged: no command line)
    QCOMPARE(f.log().back(), "PDB clean-up: wrote " + core::pathToUtf8(out));
    const core::Array entries = stepLog(f.project);
    QCOMPARE(entries.size(), std::size_t{2});
    QCOMPARE(entryText(entries[1], "step"), std::string("pdb clean"));
    QCOMPARE(entryText(entries[1], "title"), std::string("PDB clean-up of in.pdb"));
    QVERIFY(entries[1].find("command") == nullptr);
    QVERIFY(!f.window->externalChangeShown());  // a step log appended is no change of the settings
    dlg->close();
    QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    QVERIFY(!f.tools().pdbTools());
}

// The window writes what `mdbridge pdb clean` writes for the same options (the same library calls)
void ToolTests::pdbCleanWritesWhatTheCommandWrites() {
    Fixture f;
    PdbToolsWindow* dlg = f.openPdbTools();
    setText(dlg, "input", QStringLiteral("in.pdb"));
    setText(dlg, "output", QStringLiteral("window.pdb"));
    choose(dlg, "scheme", "ccd");
    setText(dlg, "l_ranges", QStringLiteral("A:2"));
    choose(dlg, "chain_mode", "list");
    setText(dlg, "chain_ids", QStringLiteral("P, Q"));
    control<QCheckBox>(dlg, "renumber")->setChecked(true);
    control<QSpinBox>(dlg, "start")->setValue(7);
    control<QCheckBox>(dlg, "per_chain")->setChecked(true);
    control<QCheckBox>(dlg, "seqres")->setChecked(true);
    setText(dlg, "fasta", QStringLiteral("window.fasta"));
    dlg->runClean();
    QVERIFY2(dlg->reportText().startsWith(QStringLiteral("Wrote window.pdb\nWrote window.fasta\n")), qPrintable(dlg->reportText()));
    QCOMPARE(runCommand({"pdb", "clean", "in.pdb", "-o", "cli.pdb", "--scheme", "ccd", "--l-ranges", "A:2", "--chain-ids",
                         "P, Q", "--renumber", "--start", "7", "--per-chain", "--seqres", "--fasta", "cli.fasta"}),
             0);
    QCOMPARE(readFile(f.dir / "window.pdb"), readFile(f.dir / "cli.pdb"));
    QCOMPARE(readFile(f.dir / "window.fasta"), readFile(f.dir / "cli.fasta"));
}

void ToolTests::builderHelixWritesWhatTheCommandWrites() {
    Fixture f;
    BuilderWindow* dlg = f.openBuilder();
    setText(dlg, "seq", QStringLiteral("GGCAUGCC"));
    choose(dlg, "chirality", "L");
    choose(dlg, "scheme", "amber");
    control<QCheckBox>(dlg, "hydrogens")->setChecked(true);
    setText(dlg, "chain1", QStringLiteral("X"));
    control<QSpinBox>(dlg, "start2")->setValue(11);
    setText(dlg, "helix_out", QStringLiteral("window.pdb"));
    click(dlg, "write");  // Watson-Crick: built at once
    QVERIFY(!dlg->busy());
    QVERIFY2(dlg->reportText().startsWith(QStringLiteral("Wrote window.pdb\n\n")), qPrintable(dlg->reportText()));
    QCOMPARE(runCommand({"build", "helix", "--seq", "GGCAUGCC", "--chirality", "L", "--scheme", "amber", "--hydrogens",
                         "--chains", "X,B", "--start", "1,11", "-o", "cli.pdb"}),
             0);
    QCOMPARE(readFile(f.dir / "window.pdb"), readFile(f.dir / "cli.pdb"));
    const core::Array entries = stepLog(f.project);
    QCOMPARE(entries.size(), std::size_t{1});
    QCOMPARE(entryText(entries[0], "title"), std::string("Ideal A-form L-RNA helix GGCAUGCC"));
    QCOMPARE(entries[0].find("citations")->get<core::Array>()->size(), std::size_t{1});  // arnott-arna
    // the chain ids as MDATK passes them (a tuple): an id typed with a comma is one long id
    setText(dlg, "chain1", QStringLiteral("A,C"));
    setText(dlg, "helix_out", QStringLiteral("comma.pdb"));
    click(dlg, "write");
    QCOMPARE(dlg->reportText(), QStringLiteral("ERROR: chain ids must be single characters"));
    QVERIFY(!fs::exists(f.dir / "comma.pdb"));
}

// MDATK tests/test_gui.py::test_builder_non_wc_helix_builds_in_the_background: the build runs on
// the window's worker and the window stays responsive; the file is what the command line writes
void ToolTests::builderNonWcHelixBuildsInTheBackground() {
    Fixture f;
    BuilderWindow* dlg = f.openBuilder();
    setText(dlg, "seq", QStringLiteral("AAAA"));
    choose(dlg, "pairing", "hh");
    const fs::path out = f.dir / "hh.pdb";
    setText(dlg, "helix_out", qs(out));
    dlg->setPage(QStringLiteral("helix"));
    dlg->runCurrent();
    QVERIFY(dlg->busy());  // the window stays responsive while it builds
    QVERIFY(!control<QPushButton>(dlg, "write")->isEnabled());
    QVERIFY(control<QProgressBar>(dlg, "progress")->isVisible());
    QVERIFY(QTest::qWaitFor([&] { return !dlg->busy(); }, 120000));
    QVERIFY(control<QPushButton>(dlg, "write")->isEnabled());
    QVERIFY(!control<QProgressBar>(dlg, "progress")->isVisible());
    const QString report = dlg->reportText();
    QVERIFY2(report.contains(QStringLiteral("Wrote")), qPrintable(report));
    QVERIFY(fs::exists(out));
    QCOMPARE(runCommand({"build", "helix", "--seq", "AAAA", "--pairing", "hh", "-o", "cli_hh.pdb"}), 0);
    QCOMPARE(readFile(out), readFile(f.dir / "cli_hh.pdb"));
    const core::Array entries = stepLog(f.project);
    QCOMPARE(entries.size(), std::size_t{1});
    QCOMPARE(entryText(entries[0], "title"), std::string("Ideal parallel tHH D-RNA helix AAAA"));
}

void ToolTests::mutateAndFitWriteWhatTheCommandsWrite() {
    Fixture f;
    BuilderWindow* dlg = f.openBuilder();
    setText(dlg, "seq", QStringLiteral("GGCAUGCC"));
    setText(dlg, "helix_out", QStringLiteral("duplex.pdb"));
    click(dlg, "write");
    dlg->setPage(QStringLiteral("mutate"));
    setText(dlg, "mut_input", QStringLiteral("duplex.pdb"));
    QCOMPARE(control<QLineEdit>(dlg, "mut_out")->text(), QStringLiteral("duplex_mut.pdb"));
    setText(dlg, "mutations", QStringLiteral("A:3G"));
    control<QCheckBox>(dlg, "wobble")->setChecked(true);
    click(dlg, "write");
    QCOMPARE(runCommand({"build", "mutate", "duplex.pdb", "A:3G", "--wobble", "-o", "cli_mut.pdb"}), 0);
    QCOMPARE(readFile(f.dir / "duplex_mut.pdb"), readFile(f.dir / "cli_mut.pdb"));

    // a map simulated from the duplex, the model fitted back into it (a local refinement)
    QCOMPARE(runCommand({"density", "simulate", "--top", "duplex.pdb", "--traj", "duplex.pdb", "--resolution", "6", "-O", "duplex.mrc"}), 0);
    dlg->setPage(QStringLiteral("fit"));
    QVERIFY(!control<QPushButton>(dlg, "check")->isEnabled());  // nothing to check on this page
    setText(dlg, "fit_input", QStringLiteral("duplex_mut.pdb"));
    QCOMPARE(control<QLineEdit>(dlg, "fit_out")->text(), QStringLiteral("duplex_mut_fit.pdb"));
    setText(dlg, "fit_map", QStringLiteral("duplex.mrc"));
    control<QDoubleSpinBox>(dlg, "resolution")->setValue(6.0);
    choose(dlg, "search", "local");
    QVERIFY(!control<QDoubleSpinBox>(dlg, "step")->isEnabled());
    click(dlg, "write");
    QVERIFY(dlg->busy());
    QVERIFY(control<QPushButton>(dlg, "cancel")->isVisible());
    QVERIFY(QTest::qWaitFor([&] { return !dlg->busy(); }, 120000));
    QVERIFY(!control<QPushButton>(dlg, "cancel")->isVisible());
    QVERIFY2(dlg->reportText().startsWith(QStringLiteral("Wrote duplex_mut_fit.pdb (CC ")), qPrintable(dlg->reportText()));
    QCOMPARE(runCommand({"build", "fit", "duplex_mut.pdb", "duplex.mrc", "--resolution", "6", "--search", "local", "-o",
                         "cli_fit.pdb"}),
             0);
    QCOMPARE(readFile(f.dir / "duplex_mut_fit.pdb"), readFile(f.dir / "cli_fit.pdb"));
    QVERIFY(QString::fromStdString(f.log().back()).startsWith(QStringLiteral("Structure builder: wrote duplex_mut_fit.pdb "
                                                                             "(density fit CC ")));
}

// TW-03 (MDATK 0.0.13): Cancel stops the fit; nothing is written, the step is logged as cancelled
void ToolTests::fitCancelIsLoggedAsCancelled() {
    Fixture f;
    BuilderWindow* dlg = f.openBuilder();
    setText(dlg, "seq", QStringLiteral("GGCAUGCCGGCAUGCC"));
    setText(dlg, "helix_out", QStringLiteral("duplex.pdb"));
    click(dlg, "write");
    QCOMPARE(runCommand({"density", "simulate", "--top", "duplex.pdb", "--traj", "duplex.pdb", "--resolution", "4", "-O", "duplex.mrc"}), 0);
    dlg->setPage(QStringLiteral("fit"));
    setText(dlg, "fit_input", QStringLiteral("duplex.pdb"));
    setText(dlg, "fit_map", QStringLiteral("duplex.mrc"));
    control<QDoubleSpinBox>(dlg, "step")->setValue(5.0);  // a long global search
    click(dlg, "write");
    QVERIFY(dlg->busy());
    QVERIFY(QTest::qWaitFor([&] { return dlg->jobs().fraction() > 0.0; }, 60000));  // the fit has begun
    click(dlg, "cancel");
    QVERIFY(QTest::qWaitFor([&] { return !dlg->busy(); }, 60000));
    QCOMPARE(dlg->reportText(), QStringLiteral("Cancelled: the density fit stopped, nothing was written."));
    QVERIFY(!fs::exists(f.dir / "duplex_fit.pdb"));
    const core::Array entries = stepLog(f.project);
    QCOMPARE(entries.size(), std::size_t{2});
    QCOMPARE(entryText(entries[1], "step"), std::string("build fit"));
    QCOMPARE(entryText(entries[1], "status"), std::string("cancelled"));
    QVERIFY(entries[1].find("files") == nullptr);
}

// TW-04 (JB-01): closing the builder while it builds asks; Cancel keeps it, Close cancels and waits
void ToolTests::closingDuringABuildAsks() {
    Fixture f;
    BuilderWindow* dlg = f.openBuilder();
    setText(dlg, "seq", QStringLiteral("GGCAUGCCGGCAUGCC"));
    setText(dlg, "helix_out", QStringLiteral("duplex.pdb"));
    click(dlg, "write");
    QCOMPARE(runCommand({"density", "simulate", "--top", "duplex.pdb", "--traj", "duplex.pdb", "--resolution", "4", "-O", "duplex.mrc"}), 0);
    dlg->setPage(QStringLiteral("fit"));
    setText(dlg, "fit_input", QStringLiteral("duplex.pdb"));
    setText(dlg, "fit_map", QStringLiteral("duplex.mrc"));
    control<QDoubleSpinBox>(dlg, "step")->setValue(5.0);
    QStringList asked;
    app::jobs::CloseAnswer answer = app::jobs::CloseAnswer::Cancel;
    dlg->jobs().setAsker([&](QWidget*, app::jobs::CloseQuestion, const QString& title) {
        asked << title;
        return answer;
    });
    click(dlg, "write");
    QVERIFY(dlg->busy());
    QVERIFY(!dlg->close());  // Cancel: the window and the job stay
    QCOMPARE(asked.size(), qsizetype{1});
    QVERIFY(asked.front().startsWith(QStringLiteral("Fitting duplex.pdb into duplex.mrc")));
    QVERIFY(dlg->isVisible() && dlg->busy());
    answer = app::jobs::CloseAnswer::Close;
    QVERIFY(dlg->close());  // cancelled and waited for
    QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    QVERIFY(!f.tools().builder());
    QVERIFY(!fs::exists(f.dir / "duplex_fit.pdb"));
}

// D1 / TW-01: the tools know only their host (milestone 7 hosts them in RNAdo)
void ToolTests::hostedWithoutTheAnalysisWindow() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    fs::copy_file(kData / "projects/default.mdatk.toml", dir / "q.mdatk.toml");
    FakeHost host;
    host.project = dir / "q.mdatk.toml";
    {
        PdbToolsWindow pdb(host);
        pdb.show();
        setText(&pdb, "input", qs(kFourResidues));
        setText(&pdb, "output", qs(dir / "out.pdb"));
        pdb.runClean();
        QCOMPARE(host.log, QStringList{QStringLiteral("PDB clean-up: wrote ") + qs(dir / "out.pdb")});
        QCOMPARE(stepLog(host.project).size(), std::size_t{1});
    }
    host.project.clear();  // an unsaved project: nothing is logged
    {
        BuilderWindow builder(host);
        builder.show();
        setText(&builder, "seq", QStringLiteral("GC"));
        setText(&builder, "helix_out", qs(dir / "gc.pdb"));
        builder.runCurrent();
        QVERIFY(fs::exists(dir / "gc.pdb"));
        QCOMPARE(host.log.back(), QStringLiteral("Structure builder: wrote ") + qs(dir / "gc.pdb"));
    }
    QCOMPARE(stepLog(dir / "q.mdatk.toml").size(), std::size_t{1});
}

// D10 as in the analysis window (WN-13): every control has an accessible name and takes the focus
// with Tab; drop-downs and spin boxes ignore the wheel unless focused
void ToolTests::everyControlIsNamedAndReachable() {
    Fixture f;
    QStringList problems;
    int checked = 0;
    for (ToolWindow* tool : std::initializer_list<ToolWindow*>{f.openPdbTools(), f.openBuilder()}) {
        for (const int page : pagesOf(tool)) {
            showPage(tool, page);
            if (auto* b = qobject_cast<BuilderWindow*>(tool)) {
                // the Cancel button shows while a fit runs: named all the same
                if (control<QPushButton>(b, "cancel")->accessibleName().isEmpty() &&
                    control<QPushButton>(b, "cancel")->text().isEmpty()) {
                    problems << QStringLiteral("builder.cancel has no name");
                }
            }
            for (QWidget* child : tool->findChildren<QWidget*>()) {
                if (!child->isVisibleTo(tool) || child->objectName().startsWith(QStringLiteral("qt_"))) continue;
                const bool focusable = (child->focusPolicy() & Qt::TabFocus) && !child->focusProxy();
                if (focusable) {
                    ++checked;
                    if (accessibleNameOf(child).isEmpty()) problems << QStringLiteral("unnamed: ") + whereIs(child);
                }
                if (isControl(child) && !qobject_cast<QScrollBar*>(child)) {
                    QWidget* owner = child;
                    while (owner->focusProxy()) owner = owner->focusProxy();
                    if (!(owner->focusPolicy() & Qt::TabFocus)) problems << QStringLiteral("no Tab focus: ") + whereIs(child);
                }
                if ((qobject_cast<QComboBox*>(child) || qobject_cast<QAbstractSpinBox*>(child)) &&
                    !app::forms::hasWheelGuard(child)) {
                    problems << QStringLiteral("the wheel changes ") + whereIs(child);
                }
            }
        }
    }
    QVERIFY2(checked > 80, qPrintable(QString::number(checked)));
    problems.removeDuplicates();
    QVERIFY2(problems.isEmpty(), qPrintable(problems.join(QLatin1Char('\n'))));
}

// Tab from a form's first control visits every enabled control of the form in the form's order
// (section by section, top to bottom), then the report and the buttons, and goes round
void ToolTests::tabOrderFollowsTheForm() {
    Fixture f;
    QStringList problems;
    for (ToolWindow* tool : std::initializer_list<ToolWindow*>{f.openPdbTools(), f.openBuilder()}) {
        for (const int page : pagesOf(tool)) {
            showPage(tool, page);
            QScrollArea* scroll = formScroll(tool, page);
            QVERIFY(scroll);
            // the form's enabled controls, top to bottom, left to right
            std::vector<QWidget*> expected;
            for (QWidget* c : scroll->widget()->findChildren<QWidget*>()) {
                if (c->isVisibleTo(tool) && c->isEnabled() && (c->focusPolicy() & Qt::TabFocus) && !c->focusProxy()) {
                    expected.push_back(c);
                }
            }
            const auto at = [&](QWidget* w) { return w->mapTo(scroll->widget(), QPoint(0, 0)); };
            std::stable_sort(expected.begin(), expected.end(), [&](QWidget* a, QWidget* b) {
                const QPoint pa = at(a), pb = at(b);
                if (std::abs(pa.y() - pb.y()) > 8) return pa.y() < pb.y();
                return pa.x() < pb.x();
            });
            QVERIFY(!expected.empty());
            tool->activateWindow();
            expected.front()->setFocus(Qt::TabFocusReason);
            QTRY_COMPARE(QApplication::focusWidget(), expected.front());
            std::vector<QWidget*> visited;
            QStringList order;
            for (int step = 0; step < 400; ++step) {
                QWidget* now = QApplication::focusWidget();
                if (!now) break;
                if (scroll->widget()->isAncestorOf(now)) visited.push_back(now);
                order << now->objectName();
                QTest::keyClick(now, Qt::Key_Tab);
                if (QApplication::focusWidget() == expected.front()) break;
            }
            if (QApplication::focusWidget() != expected.front()) {
                problems << tool->toolId() + QStringLiteral(" page %1: Tab does not go round: ").arg(page) +
                                order.join(QStringLiteral(", "));
            }
            std::vector<QWidget*> firstVisits;
            for (QWidget* w : visited) {
                if (std::find(firstVisits.begin(), firstVisits.end(), w) == firstVisits.end()) firstVisits.push_back(w);
            }
            if (firstVisits != expected) {
                QStringList want, got;
                for (QWidget* w : expected) want << w->objectName();
                for (QWidget* w : firstVisits) got << w->objectName();
                problems << tool->toolId() + QStringLiteral(" page %1: expected ").arg(page) + want.join(QStringLiteral(", ")) +
                                QStringLiteral("\n  got ") + got.join(QStringLiteral(", "));
            }
            // after the form: the report, then the buttons
            const qsizetype report = order.indexOf(tool->toolId() + QStringLiteral(".report"));
            const qsizetype close = order.indexOf(tool->toolId() + QStringLiteral(".close"));
            const qsizetype last = order.indexOf(expected.back()->objectName());
            if (!(last < report && report < close)) {
                problems << tool->toolId() + QStringLiteral(": the report and the buttons do not follow the form: ") +
                                order.join(QStringLiteral(", "));
            }
        }
    }
    QVERIFY2(problems.isEmpty(), qPrintable(problems.join(QLatin1Char('\n'))));
}

// WN-11 for the tools: at MDATK's default size (the form 520 px wide) no form is wider than its
// column (no sideways scrolling), no button or label is cut off; also at 150 % and 200 %
void ToolTests::formsFitTheirColumn() {
    Fixture f;
    QStringList overflows;
    for (ToolWindow* tool : std::initializer_list<ToolWindow*>{f.openPdbTools(), f.openBuilder()}) {
        for (const int page : pagesOf(tool)) {
            showPage(tool, page);
            QTest::qWait(50);
            QScrollArea* scroll = formScroll(tool, page);
            QWidget* viewport = scroll->viewport();
            QWidget* body = scroll->widget();
            const int room = viewport->width();
            const QString where = tool->toolId() + QStringLiteral(" page %1").arg(page);
            QVERIFY2(room >= 450, qPrintable(where + QStringLiteral(": %1 px").arg(room)));
            if (body->width() > room) {
                overflows << where + QStringLiteral(": the form is %1 px wide in %2 px").arg(body->width()).arg(room);
            }
            for (QWidget* child : body->findChildren<QWidget*>()) {
                if (!child->isVisibleTo(body)) continue;
                const auto* label = qobject_cast<QLabel*>(child);
                if ((qobject_cast<QAbstractButton*>(child) || (label && !label->wordWrap() && !label->text().isEmpty())) &&
                    (child->width() < child->minimumSizeHint().width() || child->height() < child->minimumSizeHint().height())) {
                    overflows << where + QStringLiteral(": %1 is %2 x %3 px, less than its %4 x %5 px")
                                             .arg(whereIs(child))
                                             .arg(child->width())
                                             .arg(child->height())
                                             .arg(child->minimumSizeHint().width())
                                             .arg(child->minimumSizeHint().height());
                }
                const int right = child->mapTo(viewport, QPoint(child->width(), 0)).x();
                if (right > room) overflows << where + QStringLiteral(": %1 ends at %2 px").arg(whereIs(child)).arg(right);
            }
        }
        // the button bar: nothing cut off
        for (QPushButton* b : tool->findChildren<QPushButton*>()) {
            if (b->isVisibleTo(tool) && b->width() < b->minimumSizeHint().width()) {
                overflows << whereIs(b) + QStringLiteral(" is cut off");
            }
        }
    }
    QVERIFY2(overflows.isEmpty(), qPrintable(overflows.join(QLatin1Char('\n'))));
}

void ToolTests::screenshots() {
    Fixture f;
    const QByteArray factor = qgetenv("QT_SCALE_FACTOR");
    const std::string suffix = factor.isEmpty() ? std::string() : "-" + std::to_string(static_cast<int>(factor.toDouble() * 100));
    PdbToolsWindow* pdb = f.openPdbTools();
    setText(pdb, "input", QStringLiteral("in.pdb"));
    pdb->runCheck();
    QTest::qWait(50);
    qInfo("screenshot %s", core::pathToUtf8(guitest::saveScreenshot(pdb, shotDir(), "tools-pdb_tools" + suffix)).c_str());
    BuilderWindow* builder = f.openBuilder();
    setText(builder, "seq", QStringLiteral("GGCAUGCC"));
    builder->checkCurrent();
    for (const char* page : {"helix", "mutate", "fit"}) {
        builder->setPage(QString::fromLatin1(page));
        QTest::qWait(50);
        qInfo("screenshot %s",
              core::pathToUtf8(guitest::saveScreenshot(builder, shotDir(), std::string("tools-builder-") + page + suffix))
                  .c_str());
    }
}

int main(int argc, char** argv) {
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    QApplication::setStyle(QStringLiteral("Fusion"));
    guitest::useBundledFont();
    ToolTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_tools.moc"
