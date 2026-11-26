#include "gui/actions/window_actions.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <optional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <QAbstractButton>
#include <QAction>
#include <QCoreApplication>
#include <QDesktopServices>
#include <QDir>
#include <QMetaMethod>
#include <QPushButton>
#include <QStringList>

#include "commands.hpp"
#include "gui/dialogs/image_dialog.hpp"
#include "gui/display_path.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "gui/main_window.hpp"
#include "gui/pages/input_page.hpp"
#include "gui/pages/outputs_page.hpp"
#include "gui/pages/structure_page.hpp"
#include "gui/preview_widget.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/checks.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "render_backend_qt.hpp"
#include "render_pipeline.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;

namespace {

QString qs(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }
// The output folder or one of its files with '/' separators throughout, as core::outDir resolves the
// folder (the documented Windows path rules): for the URL the desktop opens. On Windows the
// pipeline's files (core::outDir joined by std::filesystem) and a video found by listing the folder
// would otherwise end in "\\<name>". The same text on Linux and macOS.
QString generic(const fs::path& path) {
    const std::u8string text = path.generic_u8string();
    return QString::fromUtf8(reinterpret_cast<const char*>(text.data()), static_cast<qsizetype>(text.size()));
}
// ... and as the interface shows it: native separators throughout on Windows (WN-12, WA-06).
QString display(const fs::path& path) { return gui::displayPath(generic(path)); }

// What a re-render from the cache wrote.
struct Rerendered {
    std::vector<std::string> files;
    std::optional<std::string> video;
};

}  // namespace

WindowActions::WindowActions(MainWindow& window)
    : QObject(&window),
      m_window(window),
      m_backend(std::make_shared<RenderBackend>(qtRenderBackend("mdbridge"))),
      m_opener([](const QUrl& url) { return QDesktopServices::openUrl(url); }) {
    setObjectName(QStringLiteral("windowActions"));
    jobs::JobRunner& runner = *m_window.jobRunner();

    const auto wire = [&](const char* name, const QString& tip, auto slot, bool busy) {
        QAction* a = m_window.action(QString::fromLatin1(name));
        if (!a) return;
        a->setToolTip(tip);
        a->setStatusTip(tip);
        a->setEnabled(true);
        connect(a, &QAction::triggered, this, slot);
        if (busy) runner.addBusy(a);
    };
    wire("render", tr("Run cleaning, analyses and all enabled outputs (Ctrl+R)"), [this] { render(); }, true);
    wire("open_output", tr("Show the output folder in the file manager"), [this] { openOutputFolder(); }, false);
    wire("open_video", tr("Play the most recent video"), [this] { openVideo(); }, false);
    connectPages();

    connect(&window, &MainWindow::projectLoaded, this, &WindowActions::updateCentre);
    connect(window.preview(), &PreviewWidget::drawn, this, &WindowActions::updateCentre);
    connect(&window.formContext(), &forms::FormContext::documentChanged, this, &WindowActions::updateCentre);
    updateCentre();
}

WindowActions::~WindowActions() = default;

ImageDialog* WindowActions::stillDialog() const { return m_still.data(); }

void WindowActions::setUrlOpener(UrlOpener opener) {
    m_opener = opener ? std::move(opener) : UrlOpener([](const QUrl& url) { return QDesktopServices::openUrl(url); });
}

void WindowActions::connectPages() {
    jobs::JobRunner& runner = *m_window.jobRunner();
    for (const MainWindow::Page& page : m_window.pages()) {
        if (auto* input = qobject_cast<pages::InputPage*>(page.widget)) {
            connect(input, &pages::InputPage::runCleaningRequested, this, [this] { runCleaning(); });
            connect(input, &pages::InputPage::inspectSystemRequested, this, [this] { inspectSystem(); });
            runner.addBusy(input->runCleaningButton());
            runner.addBusy(input->inspectButton());
        } else if (auto* outputs = qobject_cast<pages::OutputsPage*>(page.widget)) {
            connect(outputs, &pages::OutputsPage::exportPdbsRequested, this, [this] { exportPdbs(); });
            connect(outputs, &pages::OutputsPage::rerenderRequested, this, [this] { rerenderFromCache(); });
            outputs->setJobActionsAvailable(true);
            runner.addBusy(outputs->exportPdbsButton());
            runner.addBusy(outputs->rerenderButton());
        } else if (auto* structure = qobject_cast<pages::StructurePage*>(page.widget)) {
            connect(structure, &pages::StructurePage::centreOnModelRequested, this, [this] { centreOnModel(); });
        }
        // the Layout page's Render still (milestone 5 W2-A), connected by name
        const QMetaObject* meta = page.widget->metaObject();
        const int signalIndex =
            meta->indexOfSignal(QMetaObject::normalizedSignature("renderStillRequested(int)").constData());
        if (signalIndex >= 0) {
            const int slotIndex = metaObject()->indexOfMethod(QMetaObject::normalizedSignature("onRenderStill(int)").constData());
            if (slotIndex >= 0) {
                QObject::connect(page.widget, meta->method(signalIndex), this, metaObject()->method(slotIndex));
            }
            if (auto* button = page.widget->findChild<QAbstractButton*>(QStringLiteral("layout.render_still"))) {
                runner.addBusy(button);
            }
        }
    }
}

void WindowActions::onRenderStill(int frame) { renderStill(frame); }

fs::path WindowActions::projectFile() const { return m_window.projectFile(); }

std::optional<core::Project> WindowActions::collect(bool validate, bool renderChecks, const QString& title) {
    m_window.formContext().commitPendingEdits();
    if (const QString errors = m_window.formErrorText(); !errors.isEmpty()) {
        m_window.showError(title, errors);
        return std::nullopt;
    }
    core::Project project;
    try {
        project = m_window.document().resolvedProject();
    } catch (const std::exception& e) {
        m_window.showError(title, QString::fromUtf8(e.what()));
        return std::nullopt;
    }
    if (validate) {
        std::vector<std::string> problems = core::validate(project);
        if (renderChecks) {
            for (std::string& p : render::panels::renderProblems(project)) problems.push_back(std::move(p));
        }
        if (!problems.empty()) {
            QStringList lines;
            for (const std::string& p : problems) lines.push_back(qs(p));
            m_window.showError(title, lines.join(QLatin1Char('\n')));
            return std::nullopt;
        }
    }
    return project;
}

void WindowActions::result(const QString& title, const QString& detail) {
    m_window.showInfo(title, detail);
    emit resultShown(title);
}

// ---- jobs --------------------------------------------------------------------------------------------

bool WindowActions::render() {
    std::optional<core::Project> collected = collect(true, true, tr("Please fix these parameters before rendering"));
    if (!collected) return false;
    return m_window.jobRunner()->start(
        tr("Rendering"),
        [project = std::move(*collected), file = projectFile(), backend = m_backend](const core::Reporter& reporter) {
            analysis::Pipeline analyses(project, reporter);
            RenderPipeline pipeline(analyses, *backend, file);
            RunOutputs outputs = pipeline.run(std::nullopt, std::nullopt, std::nullopt, false);
            outputs.dataset.reset();  // not needed by the window
            return outputs;
        },
        [this](RunOutputs outputs) {
            if (outputs.video) m_lastVideo = core::pathFromUtf8(*outputs.video);
            QStringList lines;
            if (outputs.video) lines.push_back(QStringLiteral("Video: %1").arg(display(core::pathFromUtf8(*outputs.video))));
            if (!outputs.charts.empty()) lines.push_back(QStringLiteral("Charts: %1 file(s)").arg(outputs.charts.size()));
            if (!outputs.pdbs.empty()) lines.push_back(QStringLiteral("PDB files: %1").arg(outputs.pdbs.size()));
            if (outputs.paramsLog) lines.push_back(QStringLiteral("Parameter log: %1").arg(display(core::pathFromUtf8(*outputs.paramsLog))));
            result(QStringLiteral("Render finished"),
                   lines.isEmpty() ? QStringLiteral("Nothing was enabled.") : lines.join(QLatin1Char('\n')));
        });
}

fs::path WindowActions::stillPath() {
    const std::string name = "MDATK_still_" + std::to_string(QCoreApplication::applicationPid()) + ".png";
    return core::pathFromUtf8(forms::toStdString(QDir::tempPath())) / name;
}

bool WindowActions::renderStill(int frame) {
    std::optional<core::Project> collected = collect(true, true, tr("Please fix these parameters first"));
    if (!collected) return false;
    const std::size_t pos = static_cast<std::size_t>(std::max(0, frame));
    return m_window.jobRunner()->start(
        tr("Rendering still"),
        [project = std::move(*collected), file = projectFile(), backend = m_backend, pos](const core::Reporter& reporter) {
            analysis::Pipeline analyses(project, reporter);
            RenderPipeline pipeline(analyses, *backend, file);
            const core::Reporter sub = reporter.sub(0.0, 0.9);
            const core::Dataset dataset = analyses.analyze(analyses.requiredAnalyses(true, false), false, &sub);
            const std::size_t last = dataset.frameCount() > 0 ? dataset.frameCount() - 1 : 0;
            return pipeline.renderPreview(dataset, std::min(pos, last), core::pathToUtf8(stillPath()));
        },
        [this](std::string written) {
            const QString path = qs(written);
            auto* dialog = new ImageDialog(&m_window, path);
            m_still = dialog;
            dialog->show();
            m_window.appendLog(QStringLiteral("Still frame written to %1").arg(path));
            emit stillShown(path);
        });
}

bool WindowActions::exportPdbs() {
    std::optional<core::Project> collected = collect(true, true, tr("Please fix these parameters first"));
    if (!collected) return false;
    return m_window.jobRunner()->start(
        tr("Exporting PDBs"),
        [project = std::move(*collected), file = projectFile(), backend = m_backend](const core::Reporter& reporter) {
            fs::create_directories(core::pathFromUtf8(core::outDir(project)));
            analysis::Pipeline analyses(project, reporter);
            RenderPipeline pipeline(analyses, *backend, file);
            return pipeline.exportPdbs();
        },
        [this](std::vector<std::string> paths) {
            QStringList lines;
            for (const std::string& p : paths) lines.push_back(display(core::pathFromUtf8(p)));
            result(QStringLiteral("PDB export finished"), lines.join(QLatin1Char('\n')));
        });
}

bool WindowActions::rerenderFromCache() {
    std::optional<core::Project> collected = collect(false, false, tr("Invalid parameters"));
    if (!collected) return false;
    return m_window.jobRunner()->start(
        tr("Re-rendering from cache"),
        [project = std::move(*collected), file = projectFile(), backend = m_backend](const core::Reporter& reporter) {
            analysis::Pipeline analyses(project, reporter);
            RenderPipeline pipeline(analyses, *backend, file);
            checkRenderSettings(project);  // as `render --from-cache` (RN-11 ... RN-15)
            fs::create_directories(core::pathFromUtf8(core::outDir(project)));
            const core::Dataset dataset = analyses.loadCachedDataset();
            Rerendered out;
            if (project.output.video.enabled) {
                const core::Reporter sub = reporter.sub(0.0, 0.85);
                out.video = pipeline.renderVideo(dataset, &sub);
                out.files.push_back(*out.video);
            }
            if (project.output.charts.enabled) {
                const core::Reporter sub = reporter.sub(0.85, 1.0);
                for (std::string& f : pipeline.renderCharts(dataset, &sub)) out.files.push_back(std::move(f));
            }
            return out;
        },
        [this](Rerendered out) {
            if (out.video) m_lastVideo = core::pathFromUtf8(*out.video);  // WA-03: the newest video
            result(QStringLiteral("Re-render finished"), QStringLiteral("%1 file(s) written").arg(out.files.size()));
        });
}

bool WindowActions::runCleaning() {
    std::optional<core::Project> collected = collect(false, false, tr("Invalid parameters"));
    if (!collected) return false;
    collected->cleaning.enabled = true;
    return m_window.jobRunner()->start(
        tr("Cleaning"),
        [project = std::move(*collected)](const core::Reporter& reporter) {
            analysis::Pipeline analyses(project, reporter);
            return analyses.resolveInputs();
        },
        [this](analysis::ResolvedInputs inputs) {
            const bool reused = inputs.cleaned && inputs.cleaned->reused;
            const QString script = inputs.cleaned ? qs(inputs.cleaned->script) : QString();
            result(reused ? QStringLiteral("Cleaned files reused") : QStringLiteral("Cleaning finished"),
                   QStringLiteral("Topology: %1\nTrajectory: %2\ncpptraj script: %3")
                       .arg(qs(inputs.topology), inputs.trajectories.empty() ? QString() : qs(inputs.trajectories.front()),
                            script));
            m_window.reloadReferenceFrame();  // WA-05: the preview shows the cleaned files
        });
}

bool WindowActions::inspectSystem() {
    std::optional<core::Project> collected = collect(false, false, tr("Invalid parameters"));
    if (!collected) return false;
    return m_window.jobRunner()->start(
        tr("Inspecting"),
        [project = std::move(*collected)](const core::Reporter& reporter) {
            const analysis::Pipeline analyses(project, reporter);
            const analysis::ResolvedInputs inputs = analyses.cleanedIfAvailable();
            std::vector<fs::path> trajectories;
            for (const std::string& t : inputs.trajectories) trajectories.push_back(core::pathFromUtf8(t));
            const core::System system =
                io::loadSystem(core::pathFromUtf8(inputs.topology), trajectories, project.inputs.dt_ps);
            for (const std::string& warning : system.warnings) reporter.log("WARNING: " + warning);
            const core::Value info = io::systemSummary(system);
            reporter.log(std::string("System (") + (inputs.topology != project.inputs.topology ? "cleaned" : "raw") +
                         " files): " + core::toJsonAscii(info, -1));
            for (const core::ResidueGroup& g : project.groups) {
                try {
                    const core::AtomGroup atoms = core::select(system, "group:" + g.name, project.groups);
                    const std::vector<std::int32_t> residues = atoms.residues();
                    if (residues.empty()) {  // WA-04 (MDATK: "ERROR index 0 is out of bounds ...")
                        reporter.log("  group " + g.name + ": 0 atoms, 0 residues");
                        continue;
                    }
                    const core::Topology& t = atoms.topology();
                    const auto label = [&t](std::int32_t r) {
                        return t.resname(r) + std::to_string(t.resids[static_cast<std::size_t>(r)]);
                    };
                    reporter.log("  group " + g.name + ": " + std::to_string(atoms.size()) + " atoms, " +
                                 std::to_string(residues.size()) + " residues (" + label(residues.front()) +
                                 "\xe2\x80\x93" + label(residues.back()) + ")");
                } catch (const std::exception& e) {  // every problem is reported, as MDATK
                    reporter.log("  group " + g.name + ": ERROR " + pythonErrorText(e));
                }
            }
        },
        [this] {
            m_window.logView()->setExpanded(true);  // WA-04: the result is in the log
            emit resultShown(QStringLiteral("Inspecting"));
        });
}

// ---- desktop -----------------------------------------------------------------------------------------

namespace {

core::Project currentProject(MainWindow& window) {
    window.formContext().commitPendingEdits();
    try {
        return window.document().resolvedProject();
    } catch (const std::exception&) {
        return window.document().project();
    }
}

}  // namespace

void WindowActions::openOutputFolder() {
    const core::Project project = currentProject(m_window);
    const fs::path folder = core::pathFromUtf8(core::outDir(project));
    std::error_code ec;
    if (!fs::is_directory(folder, ec)) {
        m_window.showError(QStringLiteral("The output directory does not exist yet."), display(folder));
        return;
    }
    if (!m_opener(QUrl::fromLocalFile(generic(folder)))) {
        m_window.showError(QStringLiteral("Could not open the folder"), tr("no file manager could be started"));
        return;
    }
    m_window.appendLog(QStringLiteral("Opened %1").arg(display(folder)));
}

fs::path WindowActions::latestVideo(const fs::path& outDir, const std::string& filename) {
    const std::string base = core::pathToUtf8(core::pathFromUtf8(filename).filename());
    const std::size_t dot = base.rfind('.');
    // os.path.splitext: a leading dot is part of the stem
    const bool hasExt = dot != std::string::npos && dot != 0 && base.find_first_not_of('.') < dot;
    const std::string stem = hasExt ? base.substr(0, dot) : base;
    const std::string ext = hasExt ? base.substr(dot) : std::string();
    fs::path best;
    fs::file_time_type bestTime{};
    std::error_code ec;
    for (fs::directory_iterator it(outDir, ec), end; !ec && it != end; it.increment(ec)) {
        const std::string name = core::pathToUtf8(it->path().filename());
        if (name.empty() || name.front() == '.') continue;  // hidden: the encoder's temporary files
        if (name.size() < stem.size() + ext.size() || name.compare(0, stem.size(), stem) != 0 ||
            name.compare(name.size() - ext.size(), ext.size(), ext) != 0) {
            continue;
        }
        std::error_code fileEc;
        if (!it->is_regular_file(fileEc)) continue;
        const fs::file_time_type time = it->last_write_time(fileEc);
        if (fileEc) continue;
        if (best.empty() || time > bestTime || (time == bestTime && it->path().filename() > best.filename())) {
            best = it->path();
            bestTime = time;
        }
    }
    return best;
}

fs::path WindowActions::videoToOpen() {
    const core::Project project = currentProject(m_window);
    const fs::path folder = core::pathFromUtf8(core::outDir(project));
    std::error_code ec;
    // WA-03: the last render's video while it is in the folder the settings name
    if (!m_lastVideo.empty() && fs::is_regular_file(m_lastVideo, ec) &&
        fs::equivalent(m_lastVideo.parent_path(), folder, ec)) {
        return m_lastVideo;
    }
    return latestVideo(folder, project.output.video.filename);
}

void WindowActions::openVideo() {
    const fs::path video = videoToOpen();
    if (video.empty()) {
        m_window.showError(QStringLiteral("No video has been rendered yet."),
                           display(core::pathFromUtf8(core::outDir(currentProject(m_window)))));
        return;
    }
    if (!m_opener(QUrl::fromLocalFile(generic(video)))) {
        m_window.showError(QStringLiteral("Could not open the video"),
                           tr("no application could be started for %1").arg(display(video)));
        return;
    }
    m_window.appendLog(QStringLiteral("Opened %1").arg(display(video)));
}

// ---- Centre on model ---------------------------------------------------------------------------------

bool WindowActions::centreAvailable() const {
    // the reference frame shown in the canvas's structure panel (W2-A), and a map
    return m_window.structurePreview() != nullptr && m_window.preview()->structureCentre().has_value() &&
           !m_window.document().project().density.map_path.empty();
}

void WindowActions::updateCentre() {
    const bool available = centreAvailable();
    for (const MainWindow::Page& page : m_window.pages()) {
        if (auto* structure = qobject_cast<pages::StructurePage*>(page.widget)) structure->setCentreOnModelAvailable(available);
    }
}

void WindowActions::centreOnModel() {
    std::optional<core::Project> collected = collect(false, false, tr("Invalid parameters"));
    if (!collected) return;
    // the structure panel of the figure the canvas shows (the example dataset on the reference frame)
    const std::optional<core::Vec3> centre =
        m_window.structurePreview() ? m_window.preview()->structureCentre() : std::nullopt;
    const std::vector<render::Vec3d>& points = m_window.preview()->structureModel();
    if (!centre) {
        m_window.showInfo(QStringLiteral("Load the real reference frame first"),
                          tr("The window loads it from the input files (see the status bar) and shows it in the "
                             "layout's structure panel; then centre the map on the model."));
        return;
    }
    std::shared_ptr<const io::DensityMap> map;
    try {
        map = io::loadMapCached(core::pathFromUtf8(collected->density.map_path));
    } catch (const std::exception& e) {
        m_window.showError(QStringLiteral("Cannot read the density map"), qs(pythonErrorText(e)));
        return;
    }
    const core::Vec3 density = render::panels::densityCentroid(*map, render::panels::finiteMapLevel(*map, collected->density.level));
    core::RigidTransform t;
    t.translation = {density[0] - (*centre)[0], density[1] - (*centre)[1], density[2] - (*centre)[2]};
    const std::vector<double> transform = t.toList();
    m_window.document().edit("Centre the map on the model", [&](core::Project& p) {
        render::panels::storeDensityPlacement(p.density, transform, points, "Centred on the model (not fitted).");
    });
}

}  // namespace mdbridge::app
