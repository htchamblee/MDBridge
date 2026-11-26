// mdbridge_gui_actions_tests: the analysis window's jobs and desktop actions (milestone 5 W2-B,
// src/app/gui/actions/window_actions.*): MDATK's test_real_frame_preview_and_render and
// test_citation_text_lists_project_resources counterparts, the enable / disable states of every
// action while a job runs, v0.0.4.1's blocking rules (D7), Export PDBs, Re-render from cache,
// Inspect system, Run cleaning, Centre on model, Open output folder / video (latest_video), Cancel
// in the middle of a render and closing the window during one (AP-14, no partial video left), the
// smoke test that makes the demo system and renders it from the window (offscreen, software
// renderer; the video parts need the FFmpeg build), and screenshots of the result dialogs (CI
// artifacts, never compared). Qt Test; messages and file dialogs answered by the harness.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QColor>
#include <QElapsedTimer>
#include <QImage>
#include <QLocale>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QTemporaryDir>
#include <QTest>
#include <QTimer>
#include <QUrl>

#include "citations.hpp"
#include "demo.hpp"
#include "gui/actions/window_actions.hpp"
#include "gui/dialogs/image_dialog.hpp"
#include "gui/dialogs/prompts.hpp"
#include "gui/dialogs/text_dialog.hpp"
#include "gui/display_path.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "gui/main_window.hpp"
#include "gui/pages/input_page.hpp"
#include "gui/pages/outputs_page.hpp"
#include "gui/pages/structure_page.hpp"
#include "gui_harness.hpp"
#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/numerics/density.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/shading.hpp"
#include "video/encoder.hpp"
#include "video/probe.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;

namespace {

constexpr int kJobWait = 300000;  // ms: a job of these tests (the demo's analyses and video)

QString qs(const fs::path& p) { return QString::fromStdString(core::pathToUtf8(p)); }
// The output folder or one of its files as the window shows it: the platform's separators
// throughout (gui::displayPath, WN-12; WindowActions' display(), WA-06), while f.outDir() / name and
// directory listings join with std::filesystem's '\\' on Windows. The same as core::pathToUtf8(p)
// on Linux and macOS.
std::string shown(const fs::path& p) { return app::gui::displayPath(p).toStdString(); }

bool videoBuild() { return render::video::videoEncodingAvailable(); }

// The files in `dir` whose names start with `prefix` and end with `suffix`.
std::vector<fs::path> filesLike(const fs::path& dir, const std::string& prefix, const std::string& suffix) {
    std::vector<fs::path> out;
    std::error_code ec;
    for (fs::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec)) {
        const std::string name = core::pathToUtf8(it->path().filename());
        if (name.size() >= prefix.size() + suffix.size() && name.rfind(prefix, 0) == 0 &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            out.push_back(it->path());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

struct Fixture {
    QTemporaryDir tmp;
    fs::path dir;
    guitest::RecordingPrompter prompter;
    QStringList nextFile;
    app::forms::FileDialogFunction previous;
    std::optional<app::MainWindow> window;
    std::vector<QUrl> opened;  // what the actions asked the desktop to open

    Fixture() {
        dir = core::pathFromUtf8(tmp.path().toStdString());
        previous = app::forms::setFileDialogFunction([this](QWidget*, const app::forms::FileDialogRequest&) {
            QStringList reply = nextFile;
            nextFile.clear();
            return reply;
        });
        window.emplace(render::Renderer::Software);
        window->setPrompter(&prompter);
        window->windowActions()->setUrlOpener([this](const QUrl& url) {
            opened.push_back(url);
            return true;
        });
        window->resize(1320, 920);
        window->show();
        if (!QTest::qWaitForWindowExposed(&*window)) qWarning("the window was not exposed");
    }
    ~Fixture() {
        if (window) {
            window->jobRunner()->stop(std::chrono::milliseconds(60000));
            window->cancelLoad();
            window->setPrompter(nullptr);
            window.reset();
        }
        app::forms::setFileDialogFunction(previous);
    }
    app::MainWindow& w() { return *window; }
    app::WindowActions& actions() { return *window->windowActions(); }
    std::vector<std::string> log() {
        window->logView()->flush();
        std::vector<std::string> out;
        for (const QString& line : window->logView()->text().split(QLatin1Char('\n'), Qt::SkipEmptyParts)) {
            out.push_back(line.toStdString());
        }
        return out;
    }
    bool logged(const std::string& line) {
        const std::vector<std::string> lines = log();
        return std::find(lines.begin(), lines.end(), line) != lines.end();
    }
    bool loggedStart(const std::string& start) {
        for (const std::string& line : log()) {
            if (line.rfind(start, 0) == 0) return true;
        }
        return false;
    }
    // MDATK's demo_dry (24 frames, no solvent) opened as MDATK's GUI fixture sets it up: output in
    // gui_out, a 320 x 180 video at 40 dpi, charts at 40 dpi of rmsd/Complex.
    fs::path openDryDemo() {
        app::DemoOptions options;
        options.frames = 24;
        options.solvent = false;
        const app::DemoFiles files = app::buildDemo(core::pathToUtf8(dir / "demo_dry"), options);
        const fs::path project = core::pathFromUtf8(files.project);
        if (!window->openPath(project, false)) return {};
        core::ProjectDocument& d = window->document();
        d.set("output.out_dir", core::pathToUtf8(dir / "gui_out"));
        d.set("output.video.width", std::int64_t{320});
        d.set("output.video.height", std::int64_t{180});
        d.set("output.video.dpi", std::int64_t{40});
        d.set("output.video.enabled", videoBuild());
        d.set("output.charts.dpi", std::int64_t{40});
        d.set("output.charts.series", core::Value(core::Array{core::Value(std::string("rmsd/Complex"))}));
        return project;
    }
    fs::path outDir() const { return dir / "gui_out"; }
    // Waits for the running job to end (events processed meanwhile); false on time-out.
    bool waitJob(int ms = kJobWait) {
        QElapsedTimer clock;
        clock.start();
        while (window->jobRunner()->running() && clock.elapsed() < ms) QTest::qWait(20);
        QApplication::processEvents();
        return !window->jobRunner()->running();
    }
    std::vector<std::string> dialogTexts(const char* kind) const {
        std::vector<std::string> out;
        for (const core::Value& d : prompter.dialogs) {
            if (d.find("kind")->stringOr() == kind) out.push_back(d.find("text")->stringOr());
        }
        return out;
    }
    std::string lastDialogText() const {
        return prompter.dialogs.empty() ? std::string() : prompter.dialogs.back().find("text")->stringOr();
    }
    std::string lastDialogDetail() const {
        return prompter.dialogs.empty() ? std::string() : prompter.dialogs.back().find("detail")->stringOr();
    }
    template <class Page>
    Page* page() {
        for (const app::MainWindow::Page& p : window->pages()) {
            if (auto* found = qobject_cast<Page*>(p.widget)) return found;
        }
        return nullptr;
    }
};

// A job that runs until it is cancelled (to look at the window while one runs).
bool startWaitingJob(app::MainWindow& w, const std::shared_ptr<std::atomic<bool>>& started) {
    return w.jobRunner()->start(
        QStringLiteral("Waiting"),
        [started](const core::Reporter& reporter) {
            started->store(true);
            for (;;) {
                reporter.checkCancel();
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
            }
        },
        [] {});
}

}  // namespace

class ActionTests : public QObject {
    Q_OBJECT

private slots:
    void initTestCase();
    void enableStates();
    void blockingRules();
    void realFramePreviewAndRender();
    void citationTextListsProjectResources();
    void aboutAndBibtexWithRunKeys();
    void exportRerenderInspect();
    void runCleaningWithoutCpptraj();
    void centreOnModel();
    void openFolderAndVideo();
    void latestVideo();
    void cancelMidRender();
    void closeDuringRender();
    void smokeDemoFromWindow();
    void resultDialogScreenshots();

private:
    QTemporaryDir m_config;
};

void ActionTests::initTestCase() {
    QVERIFY(m_config.isValid());
    qputenv("XDG_CONFIG_HOME", m_config.path().toUtf8());
    qputenv("APPDATA", m_config.path().toUtf8());
}

// Each action's enabled state, idle and while a job runs (MDATK's busy widgets and Cancel).
void ActionTests::enableStates() {
    Fixture f;
    app::MainWindow& w = f.w();
    auto* input = f.page<app::pages::InputPage>();
    auto* outputs = f.page<app::pages::OutputsPage>();
    auto* structure = f.page<app::pages::StructurePage>();
    QVERIFY(input && outputs && structure);
    const std::vector<QAbstractButton*> busyButtons{input->runCleaningButton(), input->inspectButton(),
                                                    outputs->exportPdbsButton(), outputs->rerenderButton()};
    const auto busy = [&](bool expected) {
        QCOMPARE(w.action(QStringLiteral("render"))->isEnabled(), !expected);
        QCOMPARE(w.action(QStringLiteral("create_demo"))->isEnabled(), !expected);
        QCOMPARE(w.action(QStringLiteral("cancel"))->isEnabled(), expected);
        for (QAbstractButton* b : busyButtons) QCOMPARE(b->isEnabled(), !expected);
        // opening the folder or the video never waits for a job (MDATK)
        QVERIFY(w.action(QStringLiteral("open_output"))->isEnabled());
        QVERIFY(w.action(QStringLiteral("open_video"))->isEnabled());
    };
    busy(false);
    QVERIFY(!w.action(QStringLiteral("render"))->toolTip().contains(QStringLiteral("next version")));
    // Centre on model: needs the reference frame and a map file
    QVERIFY(!structure->centreOnModelButton()->isEnabled());
    QVERIFY(!f.actions().centreAvailable());

    auto started = std::make_shared<std::atomic<bool>>(false);
    QVERIFY(startWaitingJob(w, started));
    QTRY_VERIFY(started->load());
    busy(true);
    // one job at a time: an action asked for meanwhile is refused with MDATK's message
    f.actions().rerenderFromCache();
    QCOMPARE(f.lastDialogText(), std::string("Another job is still running."));
    QVERIFY(w.jobRunner()->cancel());
    QVERIFY(f.waitJob(20000));
    busy(false);
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Cancelled"));
    QVERIFY(f.logged("Cancelled."));
}

// D7: v0.0.4.1's blocking rules (collect(validate)): typed text a field could not take blocks every
// job; Render, Render still and Export PDBs also need a project validate() accepts.
void ActionTests::blockingRules() {
    Fixture f;
    app::MainWindow& w = f.w();
    // the default settings have no input files: Render, Render still, Export PDBs refuse
    QVERIFY(!f.actions().render());
    QCOMPARE(f.lastDialogText(), std::string("Please fix these parameters before rendering"));
    QVERIFY(f.lastDialogDetail().find("topology") != std::string::npos);
    QVERIFY(!w.jobRunner()->running());
    QVERIFY(!f.actions().renderStill(0));
    QCOMPARE(f.lastDialogText(), std::string("Please fix these parameters first"));
    QVERIFY(!f.actions().exportPdbs());
    QCOMPARE(f.lastDialogText(), std::string("Please fix these parameters first"));
    QCOMPARE(f.dialogTexts("error").size(), std::size_t{3});
    // the others run without validate(): the job itself fails
    QVERIFY(f.actions().inspectSystem());
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("Inspecting failed"));
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Failed"));
    // a form error blocks them all with "Invalid parameters"
    core::Problem problem;
    problem.field = "inputs.stop";
    problem.message = "'x' is not a whole number";
    problem.location = "Input & cleaning \xe2\x80\xba Frames \xe2\x80\xba Stop";
    w.formContext().setFormError("inputs.stop", &problem);
    const std::size_t before = f.prompter.dialogs.size();
    QVERIFY(!f.actions().inspectSystem());
    QVERIFY(!f.actions().runCleaning());
    QVERIFY(!f.actions().rerenderFromCache());
    QVERIFY(!f.actions().render());
    QCOMPARE(f.prompter.dialogs.size(), before + 4);
    QCOMPARE(f.lastDialogText(), std::string("Please fix these parameters before rendering"));
    QCOMPARE(f.lastDialogDetail(), problem.text());
    QCOMPARE(f.prompter.dialogs[before].find("text")->stringOr(), std::string("Invalid parameters"));
    QVERIFY(!w.jobRunner()->running());
    w.formContext().setFormError("inputs.stop", nullptr);
}

// MDATK test_gui.py test_real_frame_preview_and_render (and its Render still).
void ActionTests::realFramePreviewAndRender() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    QVERIFY2(w.structurePreview(), qPrintable(w.loadProblem()));
    const render::StructurePreview& preview = *w.structurePreview();
    QCOMPARE(preview.trace.coords.size(), std::size_t{24});
    std::string sequence;
    for (const std::string& r : *preview.trace.resname) sequence += r;
    QCOMPARE(sequence, std::string("GCAUGGCAUCGAUCGAUGCCAUGC"));
    // the example dataset of the canvas on that frame: two chains (render.residues.chain_names)
    const core::Dataset example = render::exampleDataset(w.document().resolvedProject(), &preview);
    const std::vector<std::string> chains = render::chainNames(*example.result("trace"));
    QCOMPARE(std::set<std::string>(chains.begin(), chains.end()).size(), std::size_t{2});
    QVERIFY(w.structureSummary().contains(QStringLiteral("24 residues")));

    // Render still: one frame of the video in a dialog that fits the screen (AP-18)
    QVERIFY(f.actions().renderStill(5));
    QVERIFY(f.waitJob());
    QTRY_VERIFY(f.actions().stillDialog());
    app::ImageDialog* still = f.actions().stillDialog();
    QCOMPARE(still->image().size(), QSize(320, 180));
    QVERIFY(still->windowTitle().startsWith(QStringLiteral("MDBridge preview – MDATK_still_")));
    QCOMPARE(still->path(), qs(app::WindowActions::stillPath()));
    QVERIFY(f.logged("Still frame written to " + shown(app::WindowActions::stillPath())));
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Rendering still: done"));
    QVERIFY(QTest::qWaitForWindowExposed(still));
    guitest::saveScreenshot(still, guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR)),
                            "result-still");  // the dialog with the frame (a CI artifact)
    still->close();
    QVERIFY(f.prompter.dialogs.empty());

    // Render: one info dialog, one video, remembered for Open video
    QVERIFY(f.actions().render());
    QVERIFY(w.action(QStringLiteral("cancel"))->isEnabled());
    QVERIFY(f.waitJob());
    QCOMPARE(f.prompter.dialogs.size(), std::size_t{1});
    QCOMPARE(f.dialogTexts("info"), std::vector<std::string>{"Render finished"});
    const std::string detail = f.lastDialogDetail();
    QVERIFY(detail.find("Charts: ") != std::string::npos);
    QVERIFY(detail.find("PDB files: ") != std::string::npos);
    QVERIFY(detail.find("Parameter log: " + shown(f.outDir() / "MDATK_params_log.toml")) != std::string::npos);
    QVERIFY(f.loggedStart("Render finished: "));
    QVERIFY(f.logged("--- Rendering ---"));
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Rendering: done"));
    if (videoBuild()) {
        const std::vector<fs::path> videos = filesLike(f.outDir(), "MDATK_video_", ".mp4");
        QCOMPARE(videos.size(), std::size_t{1});
        QCOMPARE(f.actions().lastVideo(), videos.front());
        QVERIFY(detail.rfind("Video: " + shown(videos.front()), 0) == 0);
        const render::video::VideoProbe probe = render::video::probeVideo(videos.front());
        QCOMPARE(probe.width, 320);
        QCOMPARE(probe.height, 180);
        QCOMPARE(probe.packets, std::int64_t{24});
        QCOMPARE(probe.audioStreams, 0);
    }
    // the step log of the project file has the run
    w.action(QStringLiteral("project_log"))->trigger();
    QVERIFY(w.lastDialog());
    QVERIFY2(w.lastDialog()->text().contains(QStringLiteral("render  [ok]  Analyses and outputs")),
             qPrintable(w.lastDialog()->text()));
    w.lastDialog()->close();
    // the window did not take its own step log for another program's change
    w.checkDisk();
    QVERIFY(!w.externalChangeShown());

    // AP-16: About cites the settings, then what the cached results of this run record
    const app::ProjectCitations citations = app::projectCitations(w.document().resolvedProject());
    w.action(QStringLiteral("about"))->trigger();
    QCOMPARE(w.lastDialog()->text().toStdString(), app::citationText(citations));
    w.lastDialog()->close();
}

// MDATK test_gui.py test_citation_text_lists_project_resources: About's text is a non-empty text.
void ActionTests::citationTextListsProjectResources() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    w.action(QStringLiteral("about"))->trigger();
    QVERIFY(w.lastDialog());
    const QString text = w.lastDialog()->text();
    QVERIFY(!text.isEmpty());
    QVERIFY(text.startsWith(QStringLiteral("Resources this project uses (")));
    QCOMPARE(text.toStdString(), app::citationText(app::projectCitations(w.document().resolvedProject())));
    // Export BibTeX… from About writes the same keys (through withRunKeys)
    const fs::path bib = f.dir / "refs.bib";
    f.nextFile = {qs(bib)};
    w.lastDialog()->findChild<QPushButton*>(QStringLiteral("dialog.bibtex"))->click();
    QVERIFY(fs::exists(bib));
    QVERIFY(fs::file_size(bib) > 100);
}

// AP-16 (CIT-06): About and Export citations (BibTeX) add what the cached results of the last run
// record to the settings' references.
void ActionTests::aboutAndBibtexWithRunKeys() {
    Fixture f;
    app::MainWindow& w = f.w();
    w.document().set("output.out_dir", core::pathToUtf8(f.dir / "out"));
    w.document().set("cleaning.enabled", true);
    const core::Project project = w.document().resolvedProject();
    const analysis::CacheStore cache(core::cacheDir(project));
    cache.ensure();
    core::AnalysisResult result;  // a cached result recording a run-time citation (a bond guess's radii)
    result.name = "atoms";
    result.key = "0123456789abcdef";
    result.meta.emplace_back("citation_keys", core::Array{"cordero2008"});
    cache.saveResult(result, std::vector<std::int64_t>{0}, std::vector<double>{0.0});
    const core::Dataset dataset({0}, {0.0});
    cache.writeManifest(core::Table{}, core::Table{}, dataset, core::Array{}, {{"atoms", "0123456789abcdef"}});
    const app::ProjectCitations citations = app::projectCitations(project);
    QCOMPARE(citations.runOnly, std::vector<std::string>{"cordero2008"});
    w.action(QStringLiteral("about"))->trigger();
    QVERIFY(w.lastDialog());
    QVERIFY(w.lastDialog()->text().contains(QStringLiteral("\nRecorded by the cached results of the last run:\n    ")));
    QCOMPARE(w.lastDialog()->text().toStdString(), app::citationText(citations));
    const fs::path bib = f.dir / "refs.bib";
    f.nextFile = {qs(bib)};
    w.lastDialog()->findChild<QPushButton*>(QStringLiteral("dialog.bibtex"))->click();
    std::ifstream in(bib, std::ios::binary);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    QCOMPARE(text, core::citations::bibtex(citations.keys));
    QVERIFY(text.find("cordero2008") != std::string::npos);
    QVERIFY(f.logged("Wrote " + std::to_string(citations.keys.size()) + " references to " + shown(bib)));
}

void ActionTests::exportRerenderInspect() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    // Re-render before any analysis: no cache manifest
    QVERIFY(f.actions().rerenderFromCache());
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("Re-rendering from cache failed"));

    // Export PDBs now (the page's button): the demo's frames "0, last" and the middle time
    f.page<app::pages::OutputsPage>()->exportPdbsButton()->click();
    QVERIFY(w.jobRunner()->running());
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("PDB export finished"));
    const std::vector<fs::path> pdbs = filesLike(f.outDir(), "MDATK_frame_", ".pdb");
    QVERIFY(pdbs.size() >= 2);
    for (const fs::path& p : pdbs) QVERIFY(f.lastDialogDetail().find(shown(p)) != std::string::npos);

    // Re-render from the cache of a render
    w.document().set("output.pdb.enabled", false);
    QVERIFY(f.actions().render());
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("Render finished"));
    const fs::path rendered = f.actions().lastVideo();
    w.document().set("visual.title", std::string("Re-rendered"));
    QTest::qWait(1100);  // another time stamp for the new files
    f.page<app::pages::OutputsPage>()->rerenderButton()->click();
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("Re-render finished"));
    const std::size_t written = (videoBuild() ? 1u : 0u) + 1u;  // the video and the one chart
    QCOMPARE(f.lastDialogDetail(), std::to_string(written) + " file(s) written");
    if (videoBuild()) {
        QVERIFY(f.actions().lastVideo() != rendered);  // WA-03: Open video opens the re-rendered one
        QCOMPARE(f.actions().videoToOpen(), f.actions().lastVideo());
        QCOMPARE(filesLike(f.outDir(), "MDATK_video_", ".mp4").size(), std::size_t{2});
    }

    // Inspect system: the summary and the groups in the log, which unfolds (WA-04)
    w.logView()->setExpanded(false);
    const std::size_t dialogs = f.prompter.dialogs.size();
    f.page<app::pages::InputPage>()->inspectButton()->click();
    QVERIFY(f.waitJob());
    QCOMPARE(f.prompter.dialogs.size(), dialogs);  // MDATK: no dialog
    QVERIFY(w.logView()->expanded());
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Inspecting: done"));
    QVERIFY(f.loggedStart("System (raw files): {\"n_atoms\": "));
    QVERIFY(f.loggedStart("  group Aptamer: "));
    QVERIFY(f.loggedStart("  group Target: "));
    QVERIFY2(f.log().back().find(" residues (") != std::string::npos, f.log().back().c_str());
    // an empty group is reported as such
    w.document().set("groups.1.residues", std::string("900-910"));
    f.actions().inspectSystem();
    QVERIFY(f.waitJob());
    const std::vector<std::string> lines = f.log();
    const std::string target = lines.empty() ? std::string() : lines.back();
    QVERIFY2(target == "  group Target: 0 atoms, 0 residues" || target.rfind("  group Target: ERROR ", 0) == 0,
             target.c_str());
}

void ActionTests::runCleaningWithoutCpptraj() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    w.document().set("cleaning.cpptraj_path", std::string("cpptraj-not-installed"));
    f.page<app::pages::InputPage>()->runCleaningButton()->click();
    QVERIFY(w.jobRunner()->running() || !f.prompter.dialogs.empty());
    QVERIFY(f.waitJob());
    QCOMPARE(f.lastDialogText(), std::string("Cleaning failed"));
    QVERIFY(f.loggedStart("ERROR: Cleaning failed: "));
    // the project's own setting is not changed by the action (MDATK: a copy has cleaning on)
    QCOMPARE(w.document().project().cleaning.enabled, false);
}

// MDATK centre_density: the map's density centre on the structure panel's centre, no rotation; the
// placement anchored to the frame shown.
void ActionTests::centreOnModel() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    auto* structure = f.page<app::pages::StructurePage>();
    QVERIFY(!structure->centreOnModelButton()->isEnabled());  // no map yet
    // a map: one Gaussian blob around (40, -20, 15)
    io::DensityMap map;
    map.shapeXyz = {20, 20, 20};
    map.origin = {30.0, -30.0, 5.0};
    map.voxel = {1.0, 1.0, 1.0};
    map.data.resize(8000);
    for (int z = 0; z < 20; ++z) {
        for (int y = 0; y < 20; ++y) {
            for (int x = 0; x < 20; ++x) {
                const double r2 = (x - 10) * (x - 10) + (y - 10) * (y - 10) + (z - 10) * (z - 10);
                map.data[static_cast<std::size_t>((z * 20 + y) * 20 + x)] = static_cast<float>(std::exp(-r2 / 8.0));
            }
        }
    }
    const fs::path mapFile = f.dir / "blob.mrc";
    io::writeMap(mapFile, map);
    w.document().set("density.map_path", core::pathToUtf8(mapFile));
    QTRY_VERIFY(structure->centreOnModelButton()->isEnabled());
    structure->centreOnModelButton()->click();
    const core::DensityParams& d = w.document().project().density;
    QCOMPARE(d.fit_info, std::string("Centred on the model (not fitted)."));
    QCOMPARE(d.transform.size(), std::size_t{12});
    QCOMPARE(d.anchor.size(), std::size_t{24});
    QCOMPARE(d.fit_cc, 0.0);
    // rotation: identity; translation: density centre - model centre
    for (int k = 0; k < 9; ++k) QCOMPARE(d.transform[static_cast<std::size_t>(k)], k % 4 == 0 ? 1.0 : 0.0);
    const io::DensityMap read = io::readMap(mapFile);
    const core::Vec3 c = numerics::densityCentroid(read, render::panels::finiteMapLevel(read, d.level));
    QVERIFY(std::abs(c[0] - 40.0) < 1e-6 && std::abs(c[1] + 20.0) < 1e-6 && std::abs(c[2] - 15.0) < 1e-6);
    // model centre + translation = density centre: the anchor's mean is near the model centre
    QVERIFY(std::isfinite(d.transform[9]) && std::isfinite(d.transform[10]) && std::isfinite(d.transform[11]));
    QVERIFY(structure->densityStatusText().contains(QStringLiteral("Centred on the model")));
    // one undo step
    QVERIFY(w.document().canUndo());
    w.document().undo();
    QVERIFY(w.document().project().density.transform.empty());
    // a map that cannot be read
    {
        std::ofstream(f.dir / "bad.mrc") << "not a map";
    }
    w.document().set("density.map_path", core::pathToUtf8(f.dir / "bad.mrc"));
    f.actions().centreOnModel();
    QCOMPARE(f.lastDialogText(), std::string("Cannot read the density map"));
}

void ActionTests::openFolderAndVideo() {
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    // nothing written yet
    w.action(QStringLiteral("open_output"))->trigger();
    QCOMPARE(f.lastDialogText(), std::string("The output directory does not exist yet."));
    QCOMPARE(f.lastDialogDetail(), shown(f.outDir()));
    w.action(QStringLiteral("open_video"))->trigger();
    QCOMPARE(f.lastDialogText(), std::string("No video has been rendered yet."));
    QVERIFY(f.opened.empty());
    // a folder with videos of earlier runs (and the encoder's temporary file of a running one)
    fs::create_directories(f.outDir());
    const auto write = [&](const char* name, int minutesAgo) {
        const fs::path p = f.outDir() / name;
        std::ofstream(p) << name;
        fs::last_write_time(p, fs::file_time_type::clock::now() - std::chrono::minutes(minutesAgo));
        return p;
    };
    write("MDATK_video_2026-01-01_000000.mp4", 30);
    const fs::path newest = write("MDATK_video_2026-01-02_000000.mp4", 10);
    write(".MDATK_tmp_1_MDATK_video_x.mp4", 1);
    write("MDATK_video.txt", 0);
    w.action(QStringLiteral("open_output"))->trigger();
    QCOMPARE(f.opened.size(), std::size_t{1});
    QCOMPARE(f.opened.back(), QUrl::fromLocalFile(qs(f.outDir())));
    QVERIFY(f.logged("Opened " + shown(f.outDir())));
    w.action(QStringLiteral("open_video"))->trigger();
    QCOMPARE(f.opened.size(), std::size_t{2});
    QCOMPARE(f.opened.back(), QUrl::fromLocalFile(qs(newest)));
    QVERIFY(f.logged("Opened " + shown(newest)));
    // the desktop cannot open it
    w.windowActions()->setUrlOpener([](const QUrl&) { return false; });
    w.action(QStringLiteral("open_video"))->trigger();
    QCOMPARE(f.lastDialogText(), std::string("Could not open the video"));
    w.action(QStringLiteral("open_output"))->trigger();
    QCOMPARE(f.lastDialogText(), std::string("Could not open the folder"));
}

// MDATK desktop.latest_video
void ActionTests::latestVideo() {
    QTemporaryDir tmp;
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    QVERIFY(app::WindowActions::latestVideo(dir / "missing", "MDATK_video.mp4").empty());
    QVERIFY(app::WindowActions::latestVideo(dir, "MDATK_video.mp4").empty());
    const auto write = [&](const std::string& name, int minutesAgo) {
        const fs::path p = dir / name;
        std::ofstream(p) << name;
        fs::last_write_time(p, fs::file_time_type::clock::now() - std::chrono::minutes(minutesAgo));
        return p;
    };
    const fs::path plain = write("clip.mkv", 50);
    QCOMPARE(app::WindowActions::latestVideo(dir, "sub/clip.mkv"), plain);  // the base name counts
    const fs::path stamped = write("clip_2026.mkv", 40);
    write("clip_2026.mp4", 1);   // another extension
    write("other_clip.mkv", 1);  // another stem
    write(".clip_hidden.mkv", 0);
    fs::create_directories(dir / "clip_dir.mkv");  // not a file
    QCOMPARE(app::WindowActions::latestVideo(dir, "clip.mkv"), stamped);
    // equal times: the last name
    const fs::path a = write("clip_a.mkv", 5), b = write("clip_b.mkv", 5);
    fs::last_write_time(b, fs::last_write_time(a));
    QCOMPARE(app::WindowActions::latestVideo(dir, "clip.mkv"), b);
}

// Cancel in the middle of the video: no video is left, not even the encoder's temporary file.
void ActionTests::cancelMidRender() {
    if (!videoBuild()) QSKIP("built without the FFmpeg libraries");
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    w.document().set("output.video.width", std::int64_t{1280});
    w.document().set("output.video.height", std::int64_t{720});
    bool cancelled = false;
    QMetaObject::Connection c = connect(w.jobRunner(), &app::jobs::JobRunner::statusChanged, this,
                                        [&](double, const QString& text) {
                                            if (!cancelled && text.contains(QStringLiteral("Video frame"))) {
                                                cancelled = w.jobRunner()->cancel();
                                            }
                                        });
    QVERIFY(f.actions().render());
    QVERIFY(f.waitJob());
    disconnect(c);
    QVERIFY(cancelled);
    QCOMPARE(w.jobRunner()->statusText(), QStringLiteral("Cancelled"));
    QVERIFY(f.logged("Cancelled."));
    QVERIFY(f.prompter.dialogs.empty());
    QVERIFY(filesLike(f.outDir(), "", ".mp4").empty());
    QVERIFY(filesLike(f.outDir(), ".MDATK_tmp_", "").empty());
    QVERIFY(f.actions().lastVideo().empty());
}

// AP-14: closing the window during a render cancels it and waits for its clean-up.
void ActionTests::closeDuringRender() {
    if (!videoBuild()) QSKIP("built without the FFmpeg libraries");
    Fixture f;
    app::MainWindow& w = f.w();
    QVERIFY(!f.openDryDemo().empty());
    w.document().set("output.video.width", std::int64_t{1280});
    w.document().set("output.video.height", std::int64_t{720});
    QVERIFY(w.save());
    int asked = 0;
    w.jobRunner()->setAsker([&asked](QWidget*, app::jobs::CloseQuestion, const QString& title) {
        ++asked;
        return title == QStringLiteral("Rendering") ? app::jobs::CloseAnswer::Close : app::jobs::CloseAnswer::Cancel;
    });
    bool closed = false;
    QMetaObject::Connection c = connect(w.jobRunner(), &app::jobs::JobRunner::statusChanged, this,
                                        [&](double, const QString& text) {
                                            if (!closed && text.contains(QStringLiteral("Video frame"))) {
                                                closed = true;
                                                QTimer::singleShot(0, &w, [&w] { w.close(); });
                                            }
                                        });
    QVERIFY(f.actions().render());
    QTRY_VERIFY_WITH_TIMEOUT(closed && !w.jobRunner()->running(), kJobWait);
    disconnect(c);
    QCOMPARE(asked, 1);
    QVERIFY(!w.isVisible());
    QVERIFY(filesLike(f.outDir(), "", ".mp4").empty());
    QVERIFY(filesLike(f.outDir(), ".MDATK_tmp_", "").empty());
}

// The smoke test: Create demo system… and Render from the window, as a user would, with the demo's
// own settings (120 frames with solvent, 1080p video, charts, PDBs).
void ActionTests::smokeDemoFromWindow() {
    if (!videoBuild()) QSKIP("built without the FFmpeg libraries");
    Fixture f;
    app::MainWindow& w = f.w();
    const fs::path folder = f.dir / "demo";
    f.nextFile = {qs(folder)};
    w.action(QStringLiteral("create_demo"))->trigger();
    QVERIFY(f.waitJob());
    QTRY_COMPARE(w.document().path(), folder / "demo.mdatk.toml");
    QTRY_VERIFY_WITH_TIMEOUT(!w.loading(), 120000);
    QTest::keyClick(&w, Qt::Key_R, Qt::ControlModifier);  // Ctrl+R
    QVERIFY(w.jobRunner()->running());
    QVERIFY(f.waitJob());
    QCOMPARE(f.dialogTexts("error"), std::vector<std::string>{});
    QCOMPARE(f.dialogTexts("info"), std::vector<std::string>{"Render finished"});
    const fs::path out = folder / "MDATK_output";
    const std::vector<fs::path> videos = filesLike(out, "MDATK_video_", ".mp4");
    QCOMPARE(videos.size(), std::size_t{1});
    const render::video::VideoProbe probe = render::video::probeVideo(videos.front());
    QCOMPARE(probe.width, 1920);
    QCOMPARE(probe.height, 1080);
    QCOMPARE(probe.packets, std::int64_t{120});
    QCOMPARE(probe.videoStreams, 1);
    QVERIFY(probe.moovBeforeMdat);
    QVERIFY(!filesLike(out, "MDATK_frame_", ".pdb").empty());
    QVERIFY(fs::exists(out / "MDATK_params_log.toml"));
    QVERIFY(!filesLike(out, "", ".png").empty() || !filesLike(out, "", ".svg").empty());
    w.action(QStringLiteral("open_video"))->trigger();
    QCOMPARE(f.opened.back(), QUrl::fromLocalFile(qs(videos.front())));
}

// Screenshots of the result dialogs (CI artifacts, never compared): the message boxes the default
// prompter shows for them, and the still preview.
void ActionTests::resultDialogScreenshots() {
    Fixture f;
    app::MainWindow& w = f.w();
    const fs::path dir = guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR));
    const struct {
        const char* name;
        app::MessageKind kind;
        const char* text;
        std::string detail;
    } boxes[] = {
        {"render-finished", app::MessageKind::Info, "Render finished",
         "Video: /data/run/MDATK_output/MDATK_video_2026-11-26_120000.mp4\nCharts: 3 file(s)\nPDB files: 3\n"
         "Parameter log: /data/run/MDATK_output/MDATK_params_log.toml"},
        {"nothing-enabled", app::MessageKind::Info, "Render finished", "Nothing was enabled."},
        {"cleaning-finished", app::MessageKind::Info, "Cleaning finished",
         "Topology: /data/run/MDATK_cleaned/MDATK_cleaned.parm7\nTrajectory: /data/run/MDATK_cleaned/MDATK_cleaned.nc\n"
         "cpptraj script: /data/run/MDATK_cleaned/MDATK_clean.cpptraj"},
        {"pdb-export-finished", app::MessageKind::Info, "PDB export finished",
         "/data/run/MDATK_output/MDATK_frame_00000.pdb\n/data/run/MDATK_output/MDATK_frame_00119.pdb"},
        {"rerender-finished", app::MessageKind::Info, "Re-render finished", "2 file(s) written"},
        {"fix-parameters", app::MessageKind::Error, "Please fix these parameters before rendering",
         "inputs.topology is required\ninputs.trajectories needs at least one file"},
        {"rendering-failed", app::MessageKind::Error, "Rendering failed", std::string(2000, 'x')},
    };
    for (const auto& b : boxes) {
        std::unique_ptr<QMessageBox> box = app::messageBox(&w, b.kind, QString::fromUtf8(b.text), QString::fromStdString(b.detail));
        box->show();
        QVERIFY(QTest::qWaitForWindowExposed(box.get()));
        QApplication::processEvents();
        const fs::path shot = guitest::saveScreenshot(box.get(), dir, std::string("result-") + b.name);
        QVERIFY(fs::file_size(shot) > 1000);
        box->close();
    }
    // the still preview: a figure larger than the screen is scaled to fit it (AP-18)
    QImage big(4000, 2250, QImage::Format_RGB32);
    big.fill(QColor(230, 236, 245));
    const fs::path image = f.dir / "MDATK_still_1.png";
    QVERIFY(big.save(qs(image)));
    auto* still = new app::ImageDialog(&w, qs(image));
    still->show();
    QVERIFY(QTest::qWaitForWindowExposed(still));
    const QSize shown = still->shownSize();
    QVERIFY(shown.width() < 4000 && shown.height() < 2250);
    QVERIFY(std::abs(shown.width() * 2250 / 4000 - shown.height()) <= 1);
    QCOMPARE(app::ImageDialog::fitted(QSize(320, 180), QSize(1000, 800)), QSize(320, 180));
    QCOMPARE(app::ImageDialog::fitted(QSize(2000, 1000), QSize(1000, 800)), QSize(1000, 500));
    still->close();
}

int main(int argc, char** argv) {
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    QApplication::setStyle(QStringLiteral("Fusion"));
    guitest::useBundledFont();
    ActionTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_actions.moc"
