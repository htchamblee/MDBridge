// The analyze pipeline (mdatk/pipeline.py: check, resolve_inputs, universe, input_fingerprint,
// required_analyses, analyze, apply_display_overrides, load_cached_dataset; the panels' needs of
// mdatk/render/panels.py and seqpanels.py; cleaning_key and the reuse test of
// mdatk/cleaning/cpptraj.py clean()).

#include "mdbridge/analysis/pipeline.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/external.hpp"
#include "mdbridge/io/loader.hpp"
#ifdef _WIN32
#include "mdbridge/core/ntpath.hpp"
#endif

namespace mdbridge::analysis {

namespace fs = std::filesystem;
using core::Table;
using core::Value;

namespace {

// os.path.join: the folders are absolute in generic form ('/' separators) on every platform and the
// names relative, so '/' joins them.
std::string joinPath(const std::string& dir, std::string_view name) { return core::posixpath::join(dir, name); }

std::string absolute(const std::string& path) {
#ifdef _WIN32
    return core::ntpath::abspath(path, core::currentDirectory());
#else
    return core::posixpath::abspath(path, core::currentDirectory());
#endif
}

bool exists(const std::string& path) {
    std::error_code ec;
    return !path.empty() && fs::exists(core::pathFromUtf8(path), ec);
}

// Set-or-append of a table entry (a Python dict assignment).
void setEntry(Table& table, const std::string& key, Value value) {
    for (auto& [k, v] : table) {
        if (k == key) {
            v = std::move(value);
            return;
        }
    }
    table.emplace_back(key, std::move(value));
}

Table* tableEntry(Table& table, std::string_view key) {
    for (auto& [k, v] : table) {
        if (k == key) {
            return v.get<Table>();
        }
    }
    return nullptr;
}

std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

}  // namespace

// ---- cleaning -----------------------------------------------------------------------------------

Value cleaningKeyInputs(const core::CleaningParams& params, const std::string& topology,
                        const std::vector<std::string>& trajectories) {
    core::Project holder;
    holder.cleaning = params;
    Table relevant;
    const Value plain = core::toPlain(holder);
    if (const Value* section = plain.find("cleaning")) {
        for (const auto& [k, v] : *section->get<Table>()) {
            if (k != "enabled" && k != "cpptraj_path" && k != "output_dir" && k != "conda_env" && k != "conda_exe") {
                relevant.emplace_back(k, v);
            }
        }
    }
    core::Array fps;
    for (const std::string& t : trajectories) {
        fps.push_back(core::fileFingerprint(core::pathFromUtf8(t)));
    }
    const bool reference = params.fit && params.fit_reference == "file" && !params.reference_file.empty();
    return Table{{"version", 1},
                 {"topology", core::fileFingerprint(core::pathFromUtf8(topology))},
                 {"trajectories", fps},
                 {"reference", reference ? core::fileFingerprint(core::pathFromUtf8(params.reference_file)) : Value()},
                 {"params", relevant}};
}

std::string cleaningKey(const core::CleaningParams& params, const std::string& topology,
                        const std::vector<std::string>& trajectories) {
    return core::stableHash(cleaningKeyInputs(params, topology, trajectories));
}

// ---- panels -------------------------------------------------------------------------------------

std::vector<std::string> panelRequiredAnalyses(const core::Project& project) {
    std::set<std::string> need;
    const core::VisualParams& v = project.visual;
    const bool proximity = v.color_mode == "proximity" || v.glow_mode != "none";
    for (const core::PanelSpec& spec : project.layout.panels) {
        const std::string& type = spec.type;
        if (type == "structure") {
            need.insert("trace");
            if (!v.representations.empty()) {
                need.insert("atoms");
            }
            if (proximity) {
                need.insert("proximity");
            }
        } else if (type == "sequence" || type == "secondary") {
            need.insert("trace");
            if (proximity) {
                need.insert("proximity");
            }
            if (type == "secondary") {
                const Value* source = nullptr;
                for (const auto& [k, x] : spec.options) {
                    if (k == "source") {
                        source = &x;
                    }
                }
                if (!(source && source->isString() && *source->get<std::string>() == "dotbracket")) {
                    need.insert("basepairs");
                }
            }
        } else if (core::findPanelType(type) != nullptr) {
            for (const std::string& key : spec.series) {
                if (key.find('/') != std::string::npos) {
                    need.insert(key.substr(0, key.find('/')));
                }
            }
        } else {
            std::string choices;
            for (const core::PanelType& t : core::panelTypes()) {
                choices += (choices.empty() ? "" : ", ") + std::string(t.name);
            }
            throw core::ValueError("unknown panel type " + core::pythonStringRepr(type) + " (choose from " + choices +
                                   ")");
        }
    }
    return {need.begin(), need.end()};
}

// ---- Pipeline -----------------------------------------------------------------------------------

Pipeline::Pipeline(core::Project project, core::Reporter reporter)
    : m_project(std::move(project)), m_reporter(std::move(reporter)), m_cache(core::cacheDir(m_project)) {}

void Pipeline::check() const {
    const std::vector<std::string> problems = core::validate(m_project);
    if (!problems.empty()) {
        std::string text = "Invalid parameters:";
        for (const std::string& p : problems) {
            text += "\n  - " + p;
        }
        throw PipelineError(text);
    }
    for (const std::string& warning : core::layoutWarnings(m_project)) {
        m_reporter.log("WARNING: " + warning);
    }
}

const ResolvedInputs& Pipeline::resolveInputs(const core::Reporter* reporter) {
    if (m_inputs) {
        return *m_inputs;
    }
    const core::Reporter& rep = reporter ? *reporter : m_reporter;
    const core::Project& p = m_project;
    ResolvedInputs inputs;
    if (p.cleaning.enabled) {
        rep.log("Cleaning trajectory with cpptraj ...");
        // clean(): the inputs must exist, then a previous result with the same key is reused
        if (!exists(p.inputs.topology)) {
            throw core::NamedError("CpptrajError", "topology not found: " + p.inputs.topology);
        }
        for (const std::string& t : p.inputs.trajectories) {
            if (!exists(t)) {
                throw core::NamedError("CpptrajError", "trajectory not found: " + t);
            }
        }
        const std::string key = cleaningKey(p.cleaning, p.inputs.topology, p.inputs.trajectories);
        const std::string outDir = absolute(p.cleaning.output_dir.empty()
                                                ? joinPath(joinPath(core::cacheDir(p), "cleaned"), key)
                                                : p.cleaning.output_dir);
        const std::string ext = lower(core::posixpath::splitext(p.inputs.topology).second);
        const bool parm7 = ext == ".parm7" || ext == ".prmtop" || ext == ".top";
        const std::string topology = joinPath(outDir, parm7 ? kCleanedParm7 : kCleanedReference);
        const std::string trajectory = joinPath(outDir, kCleanedTrajectory);
        const std::string marker = joinPath(outDir, kCleaningMarker);
        bool reused = false;
        if (exists(marker) && exists(topology) && exists(trajectory)) {
            try {
                std::ifstream in(core::pathFromUtf8(marker), std::ios::binary);
                const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
                const Value info = core::parseJson(text, marker);
                const Value* stored = info.isTable() ? info.find("key") : nullptr;
                reused = stored && stored->isString() && *stored->get<std::string>() == key;
            } catch (const std::exception&) {
                reused = false;  // an unreadable marker counts as none
            }
        }
        if (!reused) {
            throw PipelineError("cleaning with cpptraj is not available in MDBridge yet: no cleaned trajectory "
                                "for these inputs and settings in " + outDir +
                                "; turn it off (--no-clean or cleaning.enabled = false) and give the cleaned files "
                                "as --top / --traj");
        }
        rep.rawLog("Reusing cleaned trajectory from " + outDir);
        inputs = {topology, {trajectory}, key};
    } else {
        inputs = {p.inputs.topology, p.inputs.trajectories, {}};
    }
    m_inputs = std::move(inputs);
    m_system.reset();
    rep.progress(1.0, "inputs ready");
    return *m_inputs;
}

core::System& Pipeline::system() {
    if (!m_system) {
        const ResolvedInputs& inputs = resolveInputs();
        std::vector<fs::path> trajectories;
        for (const std::string& t : inputs.trajectories) {
            trajectories.push_back(core::pathFromUtf8(t));
        }
        auto system = std::make_unique<core::System>(
            io::loadSystem(core::pathFromUtf8(inputs.topology), trajectories, m_project.inputs.dt_ps));
        for (const std::string& warning : system->warnings) {  // the readers' warnings, once
            m_reporter.log("WARNING: " + warning);
        }
        if (!system->trajectory) {
            throw core::ValueError("No trajectory loaded into Universe");
        }
        m_system = std::move(system);
    }
    return *m_system;
}

Value Pipeline::inputFingerprint(const ResolvedInputs& inputs) const {
    const core::InputParams& p = m_project.inputs;
    core::Array trajectories;
    for (const std::string& t : inputs.trajectories) {
        trajectories.push_back(core::fileFingerprint(core::pathFromUtf8(t)));
    }
    return Table{{"topology", core::fileFingerprint(core::pathFromUtf8(inputs.topology))},
                 {"trajectories", trajectories},
                 {"frames", Table{{"start", p.start}, {"stop", p.stop ? Value(*p.stop) : Value()}, {"step", p.step}}},
                 {"dt_ps", p.dt_ps ? Value(*p.dt_ps) : Value()}};
}

std::vector<std::string> Pipeline::requiredAnalyses(std::optional<bool> video, std::optional<bool> charts) const {
    const core::Project& p = m_project;
    const bool withVideo = video.value_or(p.output.video.enabled);
    const bool withCharts = charts.value_or(p.output.charts.enabled);
    std::set<std::string> names;
    for (const AnalysisDescriptor* d : registeredAnalyses()) {
        if (!d->internal() && d->enabled(p)) {
            names.insert(d->name());
        }
    }
    if (withVideo) {
        for (std::string& n : panelRequiredAnalyses(p)) {
            names.insert(std::move(n));
        }
    }
    if (withCharts) {
        for (const std::string& key : p.output.charts.series) {
            if (key.find('/') != std::string::npos) {
                names.insert(splitKey(key).first);
            }
        }
    }
    names.erase("external");
    std::vector<std::string> usable;
    for (const std::string& instance : names) {
        try {  // the affected panels / charts report "No data"
            instanceProject(p, instance);
        } catch (const std::out_of_range& e) {  // MDATK's KeyError: str() is the repr of the message
            m_reporter.log("WARNING: skipping data source " + core::pythonStringRepr(instance) + ": " +
                           core::pythonStringRepr(e.what()));
            continue;
        } catch (const core::ValueError& e) {
            m_reporter.log("WARNING: skipping data source " + core::pythonStringRepr(instance) + ": " + e.what());
            continue;
        }
        usable.push_back(instance);
    }
    std::stable_sort(usable.begin(), usable.end(), [](const std::string& a, const std::string& b) {
        const int ia = registryIndex(parseInstance(a).base);
        const int ib = registryIndex(parseInstance(b).base);
        return ia != ib ? ia < ib : a < b;
    });
    return usable;
}

core::Dataset Pipeline::analyze(const std::optional<std::vector<std::string>>& requested, bool force,
                                const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : m_reporter;
    const core::Project& p = m_project;
    const std::vector<std::string> names = requested ? *requested : requiredAnalyses();
    const core::Reporter inputsReporter = rep.sub(0.0, 0.02);
    const ResolvedInputs inputs = resolveInputs(&inputsReporter);
    const Value fp = inputFingerprint(inputs);

    // projects and keys per instance (dicts: one entry per name, first position)
    std::vector<std::string> unique;
    std::map<std::string, core::Project> projects;
    std::map<std::string, std::string> keys;
    for (const std::string& n : names) {
        if (projects.count(n) == 0) {
            projects.emplace(n, instanceProject(p, n));
            unique.push_back(n);
        }
    }
    for (const std::string& n : unique) {
        keys[n] = analysisFor(n).cacheKey(projects.at(n), fp);
    }
    m_keys.clear();
    for (const std::string& n : unique) {
        m_keys.emplace_back(n, keys[n]);
    }

    std::map<std::string, core::AnalysisResult> results;
    std::optional<std::vector<std::int64_t>> frames;
    std::vector<double> timePs;
    for (const std::string& n : names) {
        if (force) {
            break;
        }
        std::optional<CachedResult> loaded = m_cache.loadResult(n, keys[n]);
        if (loaded) {
            results[n] = std::move(loaded->result);
            frames = std::move(loaded->frames);
            timePs = std::move(loaded->timePs);
            rep.log(analysisFor(n).label() + ": loaded from cache (" + n + "-" + keys[n] + ")");
        }
    }
    std::vector<std::string> missing;
    for (const std::string& n : names) {
        if (results.count(n) == 0) {
            missing.push_back(n);
        }
    }
    if (!missing.empty()) {
        core::System& system = this->system();
        const std::vector<std::int64_t> idx =
            io::frameIndices(system.trajectory->frameCount(), p.inputs.start, p.inputs.stop, p.inputs.step);
        std::vector<std::unique_ptr<AnalysisContext>> contexts;
        std::vector<std::unique_ptr<Analysis>> analyses;
        for (const std::string& n : missing) {
            contexts.push_back(std::make_unique<AnalysisContext>(system, projects.at(n), inputs.topology, rep));
            analyses.push_back(analysisFor(n).create(*contexts.back(), n));
        }
        const core::Reporter runReporter = rep.sub(0.02, 1.0);
        RunResult run = runAnalyses(*contexts.front(), analyses, idx, &runReporter, RunOptions{m_threads});
        system.trajectory->seek(0);  // MDATK iterates trajectory[idx], which rewinds to frame 0
        for (core::AnalysisResult& res : run.results) {
            res.key = keys[res.name];
            m_cache.saveResult(res, idx, run.timePs);
            const std::string name = res.name;
            results[name] = std::move(res);
        }
        frames = idx;
        timePs = std::move(run.timePs);
    }
    if (!frames) {
        // nothing requested: still provide the frame axis
        core::System& system = this->system();
        frames = io::frameIndices(system.trajectory->frameCount(), p.inputs.start, p.inputs.stop, p.inputs.step);
        timePs.clear();
        for (const std::int64_t f : *frames) {
            timePs.push_back(system.trajectory->seek(f).time);
        }
    }
    core::Dataset dataset(*frames, timePs, {}, p.inputs.time_unit);
    for (const std::string& n : unique) {
        if (auto it = results.find(n); it != results.end()) {
            dataset.add(it->second);
        }
    }
    applyDisplayOverrides(dataset);
    if (!p.inputs.external.empty()) {
        try {
            dataset.add(io::loadExternal(p.inputs.external, dataset));
        } catch (const core::ValueError& e) {
            throw PipelineError(std::string("external data: ") + e.what());
        } catch (const fs::filesystem_error& e) {
            throw PipelineError("external data: " + core::pythonOsErrorText(e));
        }
    }
    Table frameSelection{{"start", p.inputs.start},
                         {"stop", p.inputs.stop ? Value(*p.inputs.stop) : Value()},
                         {"step", p.inputs.step}};
    core::Array trajectories(inputs.trajectories.begin(), inputs.trajectories.end());
    const Value system = Table{{"inputs", fp},
                               {"cleaning_key", inputs.cleaningKey.empty() ? Value() : Value(inputs.cleaningKey)},
                               {"topology", inputs.topology},
                               {"trajectories", trajectories}};
    core::Array groups;
    for (const core::ResidueGroup& g : p.groups) {
        groups.push_back(
            Table{{"name", g.name}, {"residues", g.residues}, {"selection", g.selection}, {"color", g.color}});
    }
    m_cache.writeManifest(system, frameSelection, dataset, groups, m_keys);
    rep.progress(1.0, "analysis done");
    return dataset;
}

void Pipeline::applyDisplayOverrides(core::Dataset& dataset) const {
    const std::optional<double>& target = m_project.analyses.energy.target_value;
    std::vector<core::AnalysisResult> changed;
    for (const core::AnalysisResult& res : dataset.results()) {
        if (parseInstance(res.name).base != "energy") {
            continue;
        }
        core::AnalysisResult copy = res;
        Table* series = tableEntry(copy.outputs(), "series");
        Table* total = series ? tableEntry(*series, "total") : nullptr;
        if (total) {
            setEntry(*total, "target", target ? Value(*target) : Value());
            changed.push_back(std::move(copy));
        }
    }
    for (core::AnalysisResult& res : changed) {
        dataset.add(std::move(res));
    }
}

core::Dataset Pipeline::loadCachedDataset() const {
    core::Dataset dataset = m_cache.loadDatasetFromManifest(m_project.inputs.time_unit);
    applyDisplayOverrides(dataset);
    return dataset;
}

}  // namespace mdbridge::analysis
