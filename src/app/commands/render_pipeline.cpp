// The output part of MDATK's pipeline (render_pipeline.hpp; mdatk/pipeline.py, MDATK's own code,
// ported).

#include "render_pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <system_error>
#include <utility>

#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/io/pdb_export.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/charts.hpp"
#include "mdbridge/render/panels/checks.hpp"
#include "mdbridge/version.hpp"
#include "density_pipeline.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;
using analysis::PipelineError;
using core::Value;

namespace {

constexpr const char* kParamsLogName = "MDATK_params_log.toml";
constexpr const char* kNoFrames = "no analysed frames to render";

[[noreturn]] void invalid(const std::vector<std::string>& problems) {
    std::string text = "Invalid parameters:";
    for (const std::string& p : problems) {
        text += "\n  - " + p;
    }
    throw PipelineError(text);
}

Value optionalText(const std::optional<std::string>& text) { return text ? Value(*text) : Value(); }

core::Array textArray(const std::vector<std::string>& items) {
    core::Array out;
    for (const std::string& s : items) out.emplace_back(s);
    return out;
}

void requireFrames(const core::Dataset& dataset) {
    if (dataset.frameCount() == 0) throw PipelineError(kNoFrames);  // RN-02
}

}  // namespace

Value RunOutputs::summary() const {
    return core::Table{{"video", optionalText(video)},
                       {"charts", textArray(charts)},
                       {"pdbs", textArray(pdbs)},
                       {"cache_manifest", optionalText(manifest)},
                       {"params_log", optionalText(paramsLog)}};
}

std::string stamped(const std::string& name, const std::string& stamp) {
    if (stamp.empty()) return name;
    const auto [root, ext] = core::posixpath::splitext(name);
    return root + "_" + stamp + ext;
}

void checkProject(const core::Project& project, const core::Reporter& reporter) {
    std::vector<std::string> problems = core::validate(project);
    for (std::string& p : render::panels::renderProblems(project)) problems.push_back(std::move(p));
    if (!problems.empty()) invalid(problems);
    for (const std::string& warning : core::layoutWarnings(project)) {
        reporter.log("WARNING: " + warning);
    }
}

void checkRenderSettings(const core::Project& project) {
    const std::vector<std::string> problems = render::panels::renderProblems(project);
    if (!problems.empty()) invalid(problems);
}

core::Table environmentInfo() {
#if defined(_WIN32)
    const char* system = "Windows";
#elif defined(__APPLE__)
    const char* system = "macOS";
#else
    const char* system = "Linux";
#endif
#if defined(__x86_64__) || defined(_M_X64)
    const char* machine = "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    const char* machine = "arm64";
#else
    const char* machine = "unknown";
#endif
    return core::Table{{"mdbridge", std::string(kVersion)},
                       {"platform", std::string(system) + " " + machine}};
}

RenderPipeline::RenderPipeline(analysis::Pipeline& analyses, const RenderBackend& backend, fs::path projectFile)
    : m_analyses(analyses), m_backend(backend), m_projectFile(std::move(projectFile)) {}

void RenderPipeline::check() const { checkProject(project(), m_analyses.reporter()); }

std::string RenderPipeline::stamp() {
    if (!project().output.timestamp_files) return {};
    if (!m_stamp) {
        const core::DateTime now = core::clockNow();
        char text[64];
        std::snprintf(text, sizeof text, "%04d-%02d-%02d_%02d%02d%02d", now.year, now.month, now.day, now.hour,
                      now.minute, now.second);
        m_stamp = text;
    }
    return *m_stamp;
}

fs::path RenderPipeline::outputPath(const std::string& name) const {
    const fs::path path = core::pathFromUtf8(name);
    return path.is_absolute() ? path : core::pathFromUtf8(core::outDir(project())) / path;
}

void RenderPipeline::startBackend() {
    if (!m_startCalled && m_backend.start) {
        m_started = m_backend.start();
        m_startCalled = true;
    }
}

const render::FontSet& RenderPipeline::fonts() {
    if (!m_backend.fonts || !m_backend.writeFigure) {
        throw std::runtime_error("figures cannot be drawn by this program (it was built without Qt); use mdbridge");
    }
    startBackend();
    return m_backend.fonts();
}

std::string RenderPipeline::renderVideo(const core::Dataset& dataset, const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : m_analyses.reporter();
    const core::Project& p = project();
    const fs::path path = outputPath(stamped(p.output.video.filename, stamp()));
    requireFrames(dataset);
    if (!p.output.video.ffmpeg_path.empty() && p.output.video.ffmpeg_path != "ffmpeg") {  // D8
        rep.log("Note: output.video.ffmpeg_path (" + p.output.video.ffmpeg_path +
                ") is not used; MDBridge encodes with its bundled FFmpeg libraries");
    }
    // ensure_density_fit: a shown map without a placement is fitted first (the first tenth of the
    // progress), MDATK's render_video
    const bool fit = core::densityActive(p) && p.density.transform.size() != 12;
    if (fit) {
        const core::Reporter first = rep.sub(0.0, 0.1);
        ensureDensityFit(m_analyses, &first);
    }
    if (!m_backend.renderVideo) {
        throw core::NamedError("VideoError", "video encoding is not available in this build");
    }
    startBackend();
    return core::pathToUtf8(m_backend.renderVideo(dataset, p, path, fit ? rep.sub(0.1, 1.0) : rep));
}

std::vector<std::string> RenderPipeline::renderCharts(const core::Dataset& dataset, const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : m_analyses.reporter();
    const core::Project& p = project();
    // the formats first (RN-11; renderCharts checks them too), then RN-02
    if (const auto problem = render::charts::chartFormatProblem(p.output.charts.formats)) {
        throw core::ValueError(*problem);
    }
    requireFrames(dataset);
    const render::FontSet& set = fonts();
    return render::charts::renderCharts(dataset, p, core::pathFromUtf8(core::outDir(p)), rep, stamp(), set,
                                        m_backend.writeFigure);
}

std::vector<std::string> RenderPipeline::exportPdbs(const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : m_analyses.reporter();
    const core::Project& p = project();
    const core::PDBParams& pp = p.output.pdb;
    core::System& system = m_analyses.system();
    const std::vector<std::int64_t> frames =
        core::framesFromSpecs(*system.trajectory, pp.frames, pp.times, p.inputs.time_unit, p.inputs.dt_ps);
    if (frames.empty()) {
        throw PipelineError("PDB export enabled but no frames or times were given");
    }
    const core::AtomGroup atoms = core::resolveExportSelection(system, pp.selection, p.groups);
    rep.log("Exporting " + std::to_string(frames.size()) + " frame(s) (" + std::to_string(atoms.size()) +
            " atoms) as PDB");
    const std::vector<fs::path> written = io::exportFrames(
        system, frames, atoms, core::pathFromUtf8(core::outDir(p)), pp.prefix,
        [&rep](double fraction, const std::string& message) { rep.progress(fraction, message); }, stamp());
    std::vector<std::string> out;
    for (const fs::path& f : written) out.push_back(core::pathToUtf8(f));
    return out;
}

std::string RenderPipeline::renderPreview(const core::Dataset& dataset, std::size_t pos, const std::string& path) {
    const core::Project& p = project();
    ensureDensityFit(m_analyses);
    requireFrames(dataset);
    if (const auto problem = render::panels::figureFileProblem(path)) throw core::ValueError(*problem);
    const std::string format = render::panels::figureFileFormat(path);
    // savefig("name") without an extension writes name.png; MDBridge prints the file it wrote
    const std::string written =
        core::posixpath::splitext(core::posixpath::basename(path)).second.empty() ? path + ".png" : path;
    const render::FontSet& set = fonts();
    render::panels::BuiltFigure built = render::panels::buildFigure(dataset, p, set);
    built.update(std::min(pos, dataset.frameCount() - 1), set);
    m_backend.writeFigure(built.figure, core::pathFromUtf8(written), format);
    return written;
}

std::string RenderPipeline::writeParamsLog(const RunOutputs& outputs) {
    const core::Project& p = project();
    const fs::path path = core::pathFromUtf8(core::outDir(p)) / kParamsLogName;
    Value inputs = core::Table{{"topology", Value()}, {"trajectories", Value()}, {"cleaning", Value()}};
    if (m_inputsResolved) {
        const analysis::ResolvedInputs& in = m_analyses.resolveInputs();
        inputs = core::Table{{"topology", in.topology},
                             {"trajectories", textArray(in.trajectories)},
                             {"cleaning", in.cleaned ? Value(core::Table{{"key", in.cleaned->key},
                                                                         {"script", in.cleaned->script},
                                                                         {"log", in.cleaned->log},
                                                                         {"reused", in.cleaned->reused}})
                                                     : Value()}};
    }
    core::Table environment = environmentInfo();
    if (m_backend.environment) {
        for (auto& entry : m_backend.environment()) environment.push_back(std::move(entry));
    }
    core::Table keys;
    for (const auto& [name, key] : m_analyses.analysisKeys()) keys.emplace_back(name, key);
    // the settings' citations, then those the cached results of this run record (CIT-06)
    const std::vector<std::string> citationKeys =
        core::citations::withRunKeys(core::citations::keysFor(core::citations::projectResources(p, {})),
                                     m_analyses.runCitationKeys());
    const core::Table run{{"project", std::string(kName)},
                          {"environment", std::move(environment)},
                          {"inputs", std::move(inputs)},
                          {"cache_dir", core::cacheDir(p)},
                          {"analysis_keys", std::move(keys)},
                          {"timestamp", optionalText(m_stamp)},
                          {"outputs", outputs.summary()},
                          {"citation_keys", textArray(citationKeys)},
                          {"citations", textArray(core::citations::textLines(citationKeys))}};
    core::writeParamsLog(path, p, run);
    return core::pathToUtf8(path);
}

RunOutputs RenderPipeline::run(std::optional<bool> video, std::optional<bool> charts, std::optional<bool> pdb,
                               bool forceAnalysis) {
    const core::Project& p = project();
    const bool v = video.value_or(p.output.video.enabled);
    const bool c = charts.value_or(p.output.charts.enabled);
    const bool d = pdb.value_or(p.output.pdb.enabled);
    check();
    fs::path logFile;  // _project_file: not a parameter log
    if (!m_projectFile.empty() && !core::isParamsLog(m_projectFile)) logFile = m_projectFile;
    RunOutputs outputs;
    const std::vector<std::string> citations =
        core::citations::keysFor(core::citations::pipelineResources(p));
    core::logged(logFile, "render", "Analyses and outputs", citations, [&](core::LogEntry& entry) {
        outputs = runSteps(v, c, d, forceAnalysis);
        entry.citations = core::citations::withRunKeys(entry.citations, m_analyses.runCitationKeys());  // CIT-06
        entry.files.clear();
        if (outputs.video) entry.files.push_back(*outputs.video);
        for (const std::string& f : outputs.charts) entry.files.push_back(f);
        for (const std::string& f : outputs.pdbs) entry.files.push_back(f);
        if (outputs.paramsLog) entry.files.push_back(*outputs.paramsLog);
        entry.summary = (outputs.dataset ? std::to_string(outputs.dataset->frameCount()) + " frames analysed"
                                         : std::string("no analyses")) +
                        "; " + std::to_string(entry.files.size()) + " file(s) written to " + core::outDir(p);
    });
    return outputs;
}

RunOutputs RenderPipeline::runSteps(bool video, bool charts, bool pdb, bool forceAnalysis) {
    const core::Project& p = project();
    const core::Reporter& rep = m_analyses.reporter();
    const std::string outDir = core::outDir(p);
    fs::create_directories(core::pathFromUtf8(outDir));
    RunOutputs outputs;
    const double start = core::clockMonotonic();
    newStamp();
    rep.log(std::string(kName) + " " + kVersion + " - output directory " + outDir);
    const core::Reporter inputsPart = rep.sub(0.0, 0.15);
    m_analyses.resolveInputs(&inputsPart);
    m_inputsResolved = true;
    const std::vector<std::string> names = m_analyses.requiredAnalyses(video, charts);
    if (!names.empty() || video || charts) {
        const core::Reporter part = rep.sub(0.15, 0.55);
        outputs.dataset = m_analyses.analyze(names, forceAnalysis, &part);
        outputs.manifest = m_analyses.cache().manifestPath();
    }
    if (video) {
        const core::Reporter part = rep.sub(0.55, 0.92);
        outputs.video = renderVideo(*outputs.dataset, &part);
    }
    if (charts) {
        const core::Reporter part = rep.sub(0.92, 0.97);
        outputs.charts = renderCharts(*outputs.dataset, &part);
    }
    if (pdb) {
        const core::Reporter part = rep.sub(0.97, 1.0);
        outputs.pdbs = exportPdbs(&part);
    }
    if (p.output.write_params_log) outputs.paramsLog = writeParamsLog(outputs);
    rep.progress(1.0, "done");
    char seconds[64];
    std::snprintf(seconds, sizeof seconds, "%.1f", core::clockMonotonic() - start);
    rep.log(std::string("Finished in ") + seconds + " s");
    return outputs;
}

}  // namespace mdbridge::app
