// MDATK's commands (mdatk/cli.py cmd_demo, cmd_init, cmd_info, cmd_list, cmd_clean, cmd_analyze,
// cmd_export_pdb, cmd_render, cmd_charts, cmd_preview, cmd_density_info and the parts of
// load_project, make_reporter and Pipeline they use; the outputs in render_pipeline.cpp), and
// MDBridge's `citations` (D14).
#include "commands.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/conda.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/io/pdb_export.hpp"
#include "mdbridge/render/panels/charts.hpp"
#include "mdbridge/render/panels/checks.hpp"
#include "chem_commands.hpp"
#include "citations.hpp"
#include "command_support.hpp"
#include "demo.hpp"
#include "render_pipeline.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;
using core::Value;

using detail::absolutePath;
using detail::addCommon;
using detail::CommonArguments;
using detail::ConsoleReporter;
using detail::lastWins;
using detail::loadProject;

struct CoreCommands::Arguments {
    CLI::App* demo = nullptr;
    CLI::App* clean = nullptr;
    CLI::App* citations = nullptr;
    CLI::App* init = nullptr;
    CLI::App* info = nullptr;
    CLI::App* list = nullptr;
    CLI::App* analyze = nullptr;
    CLI::App* exportPdb = nullptr;
    CLI::App* density = nullptr;
    CLI::App* densityInfo = nullptr;
    CLI::App* render = nullptr;
    CLI::App* charts = nullptr;
    CLI::App* preview = nullptr;
    const RenderBackend* backend = nullptr;

    std::string demoOutDir;
    std::int64_t demoFrames = 120;
    double demoDt = 10.0;
    bool demoNoSolvent = false, demoRun = false, demoQuiet = false;

    CommonArguments cleanCommon;
    bool cleanForce = false;

    CommonArguments citationsCommon;
    std::string bibtex;
    bool citationsJson = false;

    CommonArguments initCommon;
    std::string initOutput;
    bool overwrite = false;

    CommonArguments infoCommon;

    CommonArguments analyzeCommon;
    std::string only;
    bool force = false;

    CommonArguments exportCommon;
    // given (even as "") or not: CLI11 leaves an optional unset for an empty value
    std::string frameList, times, selection;
    CLI::Option* frameListOption = nullptr;
    CLI::Option* timesOption = nullptr;
    CLI::Option* selectionOption = nullptr;

    std::string map;

    CommonArguments renderCommon;
    bool video = false, noVideo = false, chartsOn = false, noCharts = false, pdb = false, noPdb = false;
    bool renderForce = false, renderFromCache = false;

    CommonArguments chartsCommon;
    std::string series, formats;
    std::int64_t dpi = 0;
    bool chartsFromCache = false;

    CommonArguments previewCommon;
    std::int64_t frame = 0;
    std::string output;
    bool previewFromCache = false;
};

namespace {


using analysis::PipelineError;

std::vector<fs::path> toPaths(const std::vector<std::string>& texts) {
    std::vector<fs::path> out;
    for (const std::string& t : texts) {
        out.push_back(core::pathFromUtf8(t));
    }
    return out;
}

// Pipeline.stamp(): the time stamp of a run's file names ("" when disabled)
std::string runStamp(const core::Project& project) {
    if (!project.output.timestamp_files) {
        return {};
    }
    const core::DateTime now = core::clockNow();
    char text[64];
    std::snprintf(text, sizeof text, "%04d-%02d-%02d_%02d%02d%02d", now.year, now.month, now.day, now.hour,
                  now.minute, now.second);
    return text;
}

int commandInit(CLI::App& command, const CommonArguments& common, const std::string& outputText, bool overwrite,
                Console& console) {
    const core::Project project = loadProject(command, common, false);  // relative paths stay relative to the file
    const fs::path output = core::pathFromUtf8(outputText);
    std::error_code ec;
    if (fs::exists(output, ec) && !overwrite) {
        console.err << outputText << " exists (use --overwrite)\n";
        return 1;
    }
    core::saveProject(output, project);
    console.out << "Wrote a project file with default settings to " << outputText << '\n';
    return 0;
}

int commandInfo(CLI::App& command, const CommonArguments& common, Console& console) {
    const core::Project project = loadProject(command, common, true);
    if (project.inputs.topology.empty()) {
        console.err << "info needs --project or --top\n";
        return 2;
    }
    const core::System system = io::loadSystem(core::pathFromUtf8(project.inputs.topology),
                                               toPaths(project.inputs.trajectories), project.inputs.dt_ps);
    for (const std::string& warning : system.warnings) {  // "WARNING: <message>" on stderr
        core::logWarning(warning);
    }
    Value info = io::systemSummary(system);
    core::Table groups;
    for (const core::ResidueGroup& g : project.groups) {
        Value entry;
        try {
            const core::AtomGroup atoms = core::select(system, "group:" + g.name, project.groups);
            const std::vector<std::int32_t> residues = atoms.residues();
            const core::Topology& t = atoms.topology();
            const auto label = [&t](std::int32_t r) {
                return t.resname(r) + std::to_string(t.resids[static_cast<std::size_t>(r)]);
            };
            entry = core::Table{{"atoms", static_cast<std::int64_t>(atoms.size())},
                                {"residues", static_cast<std::int64_t>(residues.size())},
                                {"first", label(residues.front())},
                                {"last", label(residues.back())}};
        } catch (const std::exception& error) {  // every problem is reported
            entry = core::Table{{"error", pythonErrorText(error)}};
        }
        bool replaced = false;  // a dict: a repeated name keeps its place and takes the new value
        for (auto& [name, value] : groups) {
            if (name == g.name) {
                value = entry;
                replaced = true;
            }
        }
        if (!replaced) {
            groups.emplace_back(g.name, std::move(entry));
        }
    }
    core::Array problems;
    for (const std::string& p : core::validate(project)) {
        problems.emplace_back(p);
    }
    for (const std::string& p : render::panels::renderProblems(project)) {  // RN-11 ... RN-15
        problems.emplace_back(p);
    }
    core::Table& table = *info.get<core::Table>();
    table.emplace_back("groups", std::move(groups));
    table.emplace_back("problems", std::move(problems));
    console.out << jsonForConsole(info) << '\n';
    return 0;
}

// format(text, "<width>s"): left-aligned, padded with spaces to `width` characters (code points)
std::string padded(const std::string& text, std::size_t width) {
    const std::size_t length = core::decodeUtf8Surrogateescape(text).size();
    return length >= width ? text : text + std::string(width - length, ' ');
}

int commandAnalyze(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    // _pipeline(args)
    core::Project project = loadProject(command, a.analyzeCommon, true);
    ConsoleReporter consoleReporter(console, a.analyzeCommon.quiet, a.analyzeCommon.verbose);
    analysis::Pipeline pipeline(std::move(project), consoleReporter.reporter());
    checkProject(pipeline.project(), pipeline.reporter());
    std::optional<std::vector<std::string>> names;  // args.only.split(","), each stripped
    if (!a.only.empty()) {
        names.emplace();
        std::size_t start = 0;
        while (true) {
            const std::size_t comma = a.only.find(',', start);
            const std::size_t length = comma == std::string::npos ? std::string::npos : comma - start;
            names->push_back(core::pythonStrip(std::string_view(a.only).substr(start, length)));
            if (comma == std::string::npos) {
                break;
            }
            start = comma + 1;
        }
    }
    const core::Dataset dataset = pipeline.analyze(names, a.force);
    console.out << jsonForConsole(core::Table{{"cache_manifest", pipeline.cache().manifestPath()},
                                              {"n_frames", static_cast<std::int64_t>(dataset.frameCount())},
                                              {"data_keys", dataset.describe()}})
                << '\n';
    return 0;
}

int commandExportPdb(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const CommonArguments& common = a.exportCommon;
    core::Project project = loadProject(command, common, true);
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    const core::Reporter reporter = consoleReporter.reporter();
    core::PDBParams& pdb = project.output.pdb;
    if (a.frameListOption->count() > 0) {
        pdb.frames = a.frameList;
    }
    if (a.timesOption->count() > 0) {
        pdb.times = a.times;
    }
    if (a.selectionOption->count() > 0) {
        pdb.selection = a.selection;
    }

    checkProject(project, reporter);  // Pipeline.check()
    const std::string outDir = core::outDir(project);
    fs::create_directories(core::pathFromUtf8(outDir));

    // Pipeline.export_pdbs(): universe() (resolve_inputs with cleaning, load_universe), frames, atoms, files
    analysis::Pipeline pipeline(project, reporter);
    core::System& system = pipeline.system();
    const std::vector<std::int64_t> frames = core::framesFromSpecs(
        *system.trajectory, pdb.frames, pdb.times, project.inputs.time_unit, project.inputs.dt_ps);
    if (frames.empty()) {
        throw PipelineError("PDB export enabled but no frames or times were given");
    }
    const core::AtomGroup atoms = core::resolveExportSelection(system, pdb.selection, project.groups);
    reporter.log("Exporting " + std::to_string(frames.size()) + " frame(s) (" + std::to_string(atoms.size()) +
                  " atoms) as PDB");
    const std::vector<fs::path> written = io::exportFrames(
        system, frames, atoms, core::pathFromUtf8(outDir), pdb.prefix,
        [&reporter](double fraction, const std::string& message) { reporter.progress(fraction, message); },
        runStamp(project));
    core::Array files;
    for (const fs::path& path : written) {
        files.emplace_back(core::pathToUtf8(path));
    }
    console.out << jsonForConsole(core::Table{{"pdbs", std::move(files)}}) << '\n';
    return 0;
}

// args.x.split(","), each stripped, empty ones dropped (cmd_charts)
std::vector<std::string> commaList(const std::string& text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = text.find(',', start);
        const std::string part = core::pythonStrip(
            std::string_view(text).substr(start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!part.empty()) {
            out.push_back(part);
        }
        if (comma == std::string::npos) {
            return out;
        }
        start = comma + 1;
    }
}

core::Array textArray(const std::vector<std::string>& items) {
    core::Array out;
    for (const std::string& s : items) {
        out.emplace_back(s);
    }
    return out;
}

// _dataset(pipe, args): the cached dataset (MDBridge: the render settings checked first, RN-11 ...
// RN-15), or check() and analyze()
core::Dataset commandDataset(RenderPipeline& render, bool fromCache) {
    if (fromCache) {
        checkRenderSettings(render.project());
        return render.analyses().loadCachedDataset();
    }
    render.check();
    return render.analyses().analyze();
}

int commandRender(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const CommonArguments& common = a.renderCommon;
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    analysis::Pipeline pipeline(loadProject(command, common, true), consoleReporter.reporter());
    RenderPipeline render(pipeline, *a.backend, core::pathFromUtf8(common.params));
    const core::Project& project = pipeline.project();
    if (a.renderFromCache) {
        checkRenderSettings(project);
        fs::create_directories(core::pathFromUtf8(core::outDir(project)));
        const core::Dataset dataset = pipeline.loadCachedDataset();
        core::Table out;
        if (!a.noVideo && (project.output.video.enabled || a.video)) {
            out.emplace_back("video", render.renderVideo(dataset));
        }
        if (!a.noCharts && (project.output.charts.enabled || a.chartsOn)) {
            out.emplace_back("charts", textArray(render.renderCharts(dataset)));
        }
        console.out << jsonForConsole(out) << '\n';
        return 0;
    }
    const auto choice = [](bool off, bool on) { return off ? std::optional<bool>(false) : on ? std::optional<bool>(true) : std::nullopt; };
    const RunOutputs outputs =
        render.run(choice(a.noVideo, a.video), choice(a.noCharts, a.chartsOn), choice(a.noPdb, a.pdb), a.renderForce);
    console.out << jsonForConsole(outputs.summary()) << '\n';
    return 0;
}

int commandCharts(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const CommonArguments& common = a.chartsCommon;
    core::Project project = loadProject(command, common, true);
    if (!a.series.empty()) {
        project.output.charts.series = commaList(a.series);
    }
    if (!a.formats.empty()) {
        project.output.charts.formats = commaList(a.formats);
    }
    if (a.dpi != 0) {
        project.output.charts.dpi = a.dpi;
    }
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    analysis::Pipeline pipeline(std::move(project), consoleReporter.reporter());
    RenderPipeline render(pipeline, *a.backend, core::pathFromUtf8(common.params));
    const core::Dataset dataset = commandDataset(render, a.chartsFromCache);
    console.out << jsonForConsole(core::Table{{"charts", textArray(render.renderCharts(dataset))}}) << '\n';
    return 0;
}

int commandPreview(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const CommonArguments& common = a.previewCommon;
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    analysis::Pipeline pipeline(loadProject(command, common, true), consoleReporter.reporter());
    RenderPipeline render(pipeline, *a.backend, core::pathFromUtf8(common.params));
    // RN-11: the file's format before the analyses run
    if (const auto problem = render::panels::figureFileProblem(a.output.empty() ? std::string("x.png") : a.output)) {
        throw core::ValueError(*problem);
    }
    const core::Dataset dataset = commandDataset(render, a.previewFromCache);
    if (dataset.frameCount() == 0) {  // RN-02 (MDATK: position -1, then a crash in the panels)
        throw PipelineError("no analysed frames to render");
    }
    const std::int64_t last = static_cast<std::int64_t>(dataset.frameCount()) - 1;
    const std::int64_t pos = std::min(std::max<std::int64_t>(0, a.frame), last);
    std::string out = a.output;
    if (out.empty()) {
        char name[64];
        std::snprintf(name, sizeof name, "MDATK_preview_%05lld.png", static_cast<long long>(pos));
        out = core::pathToUtf8(core::pathFromUtf8(core::outDir(pipeline.project())) / name);
    }
    fs::create_directories(core::pathFromUtf8(absolutePath(out)).parent_path());
    console.out << render.renderPreview(dataset, static_cast<std::size_t>(pos), out) << '\n';
    return 0;
}

// cmd_render with the project file alone (cmd_demo --run: _render_namespace)
int renderProject(const CLI::App& command, const std::string& projectPath, bool quiet, const RenderBackend& backend,
                  Console& console) {
    CommonArguments common;
    common.params = projectPath;
    common.quiet = quiet;
    ConsoleReporter consoleReporter(console, common.quiet, common.verbose);
    analysis::Pipeline pipeline(loadProject(command, common, true), consoleReporter.reporter());
    RenderPipeline render(pipeline, backend, core::pathFromUtf8(common.params));
    const RunOutputs outputs = render.run(std::nullopt, std::nullopt, std::nullopt, false);
    console.out << jsonForConsole(outputs.summary()) << '\n';
    return 0;
}

int commandDemo(const CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    DemoOptions options;
    options.frames = a.demoFrames;
    options.dtPs = a.demoDt;
    options.solvent = !a.demoNoSolvent;
    const DemoFiles files = buildDemo(a.demoOutDir, options);
    console.out << "Synthetic system written to " << absolutePath(a.demoOutDir) << ":\n";
    for (const auto& [key, path] : files.entries()) {
        console.out << "  " << padded(key, 10) << ' ' << path << '\n';
    }
    if (a.demoRun) {
        console.out.flush();
        return renderProject(command, files.project, a.demoQuiet, *a.backend, console);
    }
    console.out << "Next: mdbridge render -p " << files.project << "   (or: mdbridge gui " << files.project << ")\n";
    return 0;
}

int commandClean(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    core::Project project = loadProject(command, a.cleanCommon, true);
    project.cleaning.enabled = true;
    ConsoleReporter consoleReporter(console, a.cleanCommon.quiet, a.cleanCommon.verbose);
    analysis::Pipeline pipeline(std::move(project), consoleReporter.reporter());
    const analysis::ResolvedInputs& res = pipeline.resolveInputs(nullptr, a.cleanForce);
    const std::optional<analysis::CleanResult>& c = res.cleaned;
    console.out << jsonForConsole(core::Table{{"topology", res.topology},
                                              {"trajectory", textArray(res.trajectories)},
                                              {"script", c ? Value(c->script) : Value()},
                                              {"log", c ? Value(c->log) : Value()},
                                              {"reused", c ? Value(c->reused) : Value()}})
                << '\n';
    return 0;
}

// MDBridge's citations (D14; tests/parity/citations_command.py)
int commandCitations(CLI::App& command, const CoreCommands::Arguments& a, Console& console) {
    const core::Project project = loadProject(command, a.citationsCommon, true);
    const ProjectCitations citations = projectCitations(project);
    if (!a.bibtex.empty()) {
        if (citations.keys.empty()) {
            console.out << "Nothing to cite yet: The project uses no external resource with a reference.\n";
            return 0;
        }
        std::ofstream out(core::pathFromUtf8(a.bibtex), std::ios::binary | std::ios::trunc);
        if (!out) {
            throw fs::filesystem_error("cannot write the BibTeX file", core::pathFromUtf8(a.bibtex),
                                       std::error_code(errno, std::generic_category()));
        }
        out << core::citations::bibtex(citations.keys);
        out.close();
        if (!out) {
            throw fs::filesystem_error("cannot write the BibTeX file", core::pathFromUtf8(a.bibtex),
                                       std::make_error_code(std::errc::io_error));
        }
    }
    if (a.citationsJson) {
        console.out << core::toJson(citationJson(citations), 2) << '\n';
    } else if (!a.bibtex.empty()) {
        console.out << "Wrote " << citations.keys.size() << " references to " << a.bibtex << '\n';
    } else {
        console.out << citationText(citations);
    }
    return 0;
}

int commandList(Console& console) {
    printList(console.out);
    return 0;
}

int commandDensityInfo(const std::string& map, Console& console) {
    console.out << jsonForConsole(io::densityInfo(core::pathFromUtf8(map))) << '\n';
    return 0;
}

}  // namespace

// ---- the command line ------------------------------------------------------------------------

Console standardConsole() {
#ifdef _WIN32
    const bool tty = _isatty(_fileno(stderr)) != 0;
#else
    const bool tty = isatty(fileno(stderr)) != 0;
#endif
    return Console{std::cout, std::cerr, tty};
}

CoreCommands::CoreCommands(CLI::App& app, const RenderBackend* backend, const NetworkBackend* network)
    : m_args(std::make_unique<Arguments>()) {
    Arguments& a = *m_args;
    static const RenderBackend kNoBackend;
    a.backend = backend ? backend : &kNoBackend;

    a.demo = app.add_subcommand("demo", "write a synthetic RNA aptamer/target system and project file");
    a.demo->add_option("out_dir", a.demoOutDir, "folder for demo.parm7, demo.nc and demo.mdatk.toml")->required();
    lastWins(a.demo->add_option("--n-frames", a.demoFrames, "number of frames"))->capture_default_str();
    lastWins(a.demo->add_option("--dt", a.demoDt, "ps between frames"))->capture_default_str();
    a.demo->add_flag("--no-solvent", a.demoNoSolvent, "omit water / Na+ (no cleaning needed)");
    a.demo->add_flag("--run", a.demoRun, "render the demo immediately");
    a.demo->add_flag("--quiet,-q", a.demoQuiet, "only print the result");

    a.init = app.add_subcommand("init", "write a project file with default settings");
    a.init->alias("init-params");
    a.init->add_option("output", a.initOutput, "project file to write (.mdatk.toml)")->required();
    a.init->add_flag("--overwrite", a.overwrite, "replace an existing file");
    addCommon(*a.init, a.initCommon);

    a.info = app.add_subcommand("info", "summarise a system and check residue groups");
    addCommon(*a.info, a.infoCommon);

    a.list = app.add_subcommand("list", "list analyses, panel types and layout presets");

    a.clean = app.add_subcommand("clean", "run (or reuse) cpptraj stripping/imaging/fitting");
    addCommon(*a.clean, a.cleanCommon);
    a.clean->add_flag("--force", a.cleanForce, "re-run cpptraj even if cached");

    a.analyze = app.add_subcommand("analyze", "compute (or load cached) analyses");
    addCommon(*a.analyze, a.analyzeCommon);
    lastWins(a.analyze->add_option("--only", a.only, "comma separated analysis names"));
    a.analyze->add_flag("--force", a.force, "ignore cached results");

    a.exportPdb = app.add_subcommand("export-pdb", "write selected frames as PDB files");
    addCommon(*a.exportPdb, a.exportCommon);
    a.frameListOption =
        lastWins(a.exportPdb->add_option("--frame-list", a.frameList, "frames, e.g. \"0, 10, 20-50:10, last\""));
    a.timesOption = lastWins(a.exportPdb->add_option("--times", a.times, "times, e.g. \"1.5ns, 2000ps\""));
    a.selectionOption = lastWins(
        a.exportPdb->add_option("--selection", a.selection, "group names (comma separated) or a selection"));

    a.render = app.add_subcommand("render", "full run: clean, analyse, video, charts, PDBs");
    addCommon(*a.render, a.renderCommon);
    a.render->add_flag("--video", a.video, "force video output on");
    a.render->add_flag("--no-video", a.noVideo, "skip video output");
    a.render->add_flag("--charts", a.chartsOn, "force charts output on");
    a.render->add_flag("--no-charts", a.noCharts, "skip charts output");
    a.render->add_flag("--pdb", a.pdb, "force pdb output on");
    a.render->add_flag("--no-pdb", a.noPdb, "skip pdb output");
    a.render->add_flag("--force", a.renderForce, "recompute analyses");
    a.render->add_flag("--from-cache", a.renderFromCache,
                       "re-render from the cache manifest without touching trajectories");

    a.charts = app.add_subcommand("charts", "write standalone time charts");
    addCommon(*a.charts, a.chartsCommon);
    lastWins(a.charts->add_option("--series", a.series, "comma separated data keys (default: all)"));
    lastWins(a.charts->add_option("--formats", a.formats, "comma separated: png,svg,pdf"));
    lastWins(a.charts->add_option("--dpi", a.dpi, "raster resolution"));
    a.charts->add_flag("--from-cache", a.chartsFromCache, "use the analyses of the cache manifest");

    a.preview = app.add_subcommand("preview", "render one composite frame as an image");
    addCommon(*a.preview, a.previewCommon);
    lastWins(a.preview->add_option("--frame", a.frame, "position in the analysed frames"))->capture_default_str();
    lastWins(a.preview->add_option("--output", a.output, "image path (.png/.svg/.pdf)"));
    a.preview->add_flag("--from-cache", a.previewFromCache, "use the analyses of the cache manifest");

    a.citations = app.add_subcommand("citations", "list the references a project's settings and results cite");
    addCommon(*a.citations, a.citationsCommon);
    lastWins(a.citations->add_option("--bibtex", a.bibtex, "write the references as BibTeX to FILE"))->type_name("FILE");
    a.citations->add_flag("--json", a.citationsJson, "print the keys and resources as JSON");

    a.density = app.add_subcommand("density", "density maps (MRC / CCP4): header info, rigid-body fit of the "
                                              "model, simulated maps");
    a.density->require_subcommand(0, 1);
    a.densityInfo = a.density->add_subcommand("info", "print a map's grid, voxel size, origin and statistics");
    a.densityInfo->add_option("map", a.map, "MRC / CCP4 map (.mrc, .map, .ccp4; gzip-compressed too)")->required();

    // density fit / simulate, pdb, build, dict (milestone 6)
    m_chem = std::make_unique<detail::ChemCommands>(app, *a.density, network);
}

CoreCommands::~CoreCommands() = default;

std::optional<int> CoreCommands::run(Console& console) {
    Arguments& a = *m_args;
    const bool any = a.demo->parsed() || a.clean->parsed() || a.citations->parsed() || a.init->parsed() || a.info->parsed() || a.list->parsed() || a.analyze->parsed() ||
                     a.exportPdb->parsed() || a.density->parsed() || a.render->parsed() || a.charts->parsed() ||
                     a.preview->parsed() || m_chem->parsed();
    if (!any) {
        return std::nullopt;
    }
    if (a.density->parsed() && !a.densityInfo->parsed() && !m_chem->densityParsed()) {
        console.out << a.density->help();
        return 2;
    }
    // logging.basicConfig(level=WARNING, format="%(levelname)s: %(message)s")
    const core::LogSink previous = core::setLogSink([&console](core::LogLevel level, std::string_view message) {
        if (level == core::LogLevel::Warning || level == core::LogLevel::Error) {
            console.err << core::logLevelName(level) << ": " << message << '\n';
            console.err.flush();
        }
    });
    struct RestoreSink {
        const core::LogSink& sink;
        ~RestoreSink() { core::setLogSink(sink); }
    } restore{previous};
    try {
        if (a.demo->parsed()) {
            return commandDemo(*a.demo, a, console);
        }
        if (a.clean->parsed()) {
            return commandClean(*a.clean, a, console);
        }
        if (a.citations->parsed()) {
            return commandCitations(*a.citations, a, console);
        }
        if (a.init->parsed()) {
            return commandInit(*a.init, a.initCommon, a.initOutput, a.overwrite, console);
        }
        if (a.info->parsed()) {
            return commandInfo(*a.info, a.infoCommon, console);
        }
        if (a.list->parsed()) {
            return commandList(console);
        }
        if (a.analyze->parsed()) {
            return commandAnalyze(*a.analyze, a, console);
        }
        if (a.exportPdb->parsed()) {
            return commandExportPdb(*a.exportPdb, a, console);
        }
        if (a.render->parsed()) {
            return commandRender(*a.render, a, console);
        }
        if (a.charts->parsed()) {
            return commandCharts(*a.charts, a, console);
        }
        if (a.preview->parsed()) {
            return commandPreview(*a.preview, a, console);
        }
        if (m_chem->parsed()) {
            return m_chem->run(console, m_arguments);
        }
        return commandDensityInfo(a.map, console);
    } catch (const core::Cancelled&) {
        console.err << "cancelled\n";
        return 130;
    } catch (const std::exception& error) {
        console.out.flush();
        console.err << "error: " << pythonErrorText(error) << '\n';
        return 1;
    }
}

int runCoreCommandLine(int argc, const char* const* argv, Console& console, const RenderBackend* backend,
                       const NetworkBackend* network) {
    CLI::App app("MDBridge: the MDATK commands of the core and io libraries (no interface)",
                 argc > 0 ? core::pathBasename(argv[0]) : "mdbridge");
    app.require_subcommand(0, 1);
    CoreCommands commands(app, backend, network);
    commands.setArguments(std::vector<std::string>(argv + std::min(argc, 1), argv + argc));
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        const int code = app.exit(error, console.out, console.err);
        return code == 0 ? 0 : 2;
    }
    if (std::optional<int> code = commands.run(console)) {
        return *code;
    }
    console.out << app.help();
    return 2;
}

// ---- list -------------------------------------------------------------------------------------

void printList(std::ostream& out) {
    out << "Analyses:\n";
    for (const analysis::AnalysisDescriptor* d : analysis::registeredAnalyses()) {
        out << "  " << padded(d->name(), 10) << ' ' << d->label()
            << (d->internal() ? " (internal, used by the structure panel)" : "") << '\n';
    }
    out << "\nConda environments with cpptraj:\n";
    const std::vector<core::conda::CondaEnv> envs = core::conda::findCpptrajEnvs();
    for (const core::conda::CondaEnv& env : envs) {
        out << "  " << padded(env.name, 14) << ' ' << env.prefix << '\n';
    }
    if (envs.empty()) {
        out << "  (none found)\n";
    }
    out << "\nPanel types:\n";
    for (const core::PanelType& type : core::panelTypes()) {
        out << "  " << padded(std::string(type.name), 12) << ' ' << type.label << " - " << type.description << '\n';
    }
    out << "\nLayout presets:\n";
    for (const core::LayoutPreset& preset : core::layoutPresets()) {
        out << "  " << padded(std::string(preset.name), 14) << ' ' << preset.label << '\n';
    }
    out << "\nData keys for panels/charts look like '<analysis>/<output>', e.g. rmsd/Complex, contacts/atom_pairs,\n"
           "contacts/per_residue_a (matrix), contacts/map (contact map). Add residue groups after '@' to\n"
           "run an analysis on other groups: contacts@Target+Mg/atom_pairs, rmsd@Aptamer/Backbone.\n"
           "Run 'mdbridge analyze' to list the keys of a system.\n";
    out.flush();
}

// ---- Python's texts --------------------------------------------------------------------------

std::string jsonForConsole(const Value& value) { return core::toJsonAscii(value, 2); }

std::string pythonErrorText(const std::exception& error) {
    if (dynamic_cast<const fs::filesystem_error*>(&error) != nullptr) {
        // io's FileNotFoundError("topology not found: <path>") carries its message alone
        const std::string what = error.what();
        if (what.starts_with("topology not found: ") || what.starts_with("trajectory not found: ")) {
            return what;
        }
    }
    if (dynamic_cast<const std::out_of_range*>(&error) != nullptr) {
        // MDATK prints str(KeyError), the repr of the message ("'...'" / "\"...\""): the message alone
        return error.what();
    }
    const std::string summary = core::exceptionSummary(error);  // "<Python class>: <str(exc)>"
    const std::size_t colon = summary.find(": ");
    return colon == std::string::npos ? std::string() : summary.substr(colon + 2);
}

}  // namespace mdbridge::app
