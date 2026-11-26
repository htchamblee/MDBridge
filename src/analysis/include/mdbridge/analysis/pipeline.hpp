#pragma once
// The analyze part of MDATK's end-to-end pipeline (mdatk/pipeline.py Pipeline): checking a
// project, resolving its input files, the analyses a run needs (enabled analyses, the panels'
// data and the charted series), running the missing ones in one pass and reusing the rest from the
// cache (cache.hpp), the dataset with display overrides and external data, and the manifest.
//
// With cleaning.enabled the inputs are first cleaned with cpptraj (cleaning.hpp), or the cleaned
// files of an earlier run with the same inputs and settings (MDATK's or MDBridge's) are reused.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/analysis/cleaning.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::analysis {

// pipeline.PipelineError (a RuntimeError): what the CLI reports as "error: <message>".
class PipelineError : public core::NamedError {
public:
    explicit PipelineError(const std::string& message) : core::NamedError("PipelineError", message) {}
};

// The files MDATK's cleaning writes (mdatk/cleaning/cpptraj.py).
inline constexpr std::string_view kCleanedParm7 = "MDATK_cleaned.parm7";
inline constexpr std::string_view kCleanedTrajectory = "MDATK_cleaned.nc";
inline constexpr std::string_view kCleanedReference = "MDATK_cleaned_ref.pdb";
inline constexpr std::string_view kCleaningMarker = "MDATK_clean.json";

// What cleaning_key(params, topology, trajectories) hashes: {"version": 1, "topology":
// fileFingerprint, "trajectories": [...], "reference": fingerprint of cleaning.reference_file when
// fit with fit_reference "file" (else null), "params": the cleaning settings without enabled,
// cpptraj_path, output_dir, conda_env and conda_exe}. Throws std::filesystem::filesystem_error for
// a missing file.
core::Value cleaningKeyInputs(const core::CleaningParams& params, const std::string& topology,
                              const std::vector<std::string>& trajectories);
// cleaning_key: stableHash(cleaningKeyInputs(...)), 16 hex digits.
std::string cleaningKey(const core::CleaningParams& params, const std::string& topology,
                        const std::vector<std::string>& trajectories);

// render.panels.required_analyses(project): the analysis instances the layout's panels read.
// structure: trace (+ atoms with representations, + proximity when visual.color_mode is
// "proximity" or visual.glow_mode is not "none"); sequence: trace (+ proximity); secondary: trace
// (+ proximity) + basepairs unless option source is "dotbracket"; the other types: the instance
// part of their series keys that contain "/". Sorted. Throws core::ValueError("unknown panel type
// 'x' (choose from structure, ...)") for a type core::panelTypes() lacks.
std::vector<std::string> panelRequiredAnalyses(const core::Project& project);

struct ResolvedInputs {
    std::string topology;
    std::vector<std::string> trajectories;
    std::string cleaningKey;  // the cleaned files' key; "" without cleaning
    // clean()'s result (script, log, reused) when resolveInputs cleaned or reused cleaned files
    std::optional<CleanResult> cleaned;
};

class Pipeline {
public:
    // The project as loaded (paths resolved); the reporter receives log lines and progress.
    explicit Pipeline(core::Project project, core::Reporter reporter = core::Reporter());

    const core::Project& project() const { return m_project; }
    const core::Reporter& reporter() const { return m_reporter; }
    const CacheStore& cache() const { return m_cache; }
    // Worker threads of analyze()'s run (RunOptions::threads: 0, the default, for
    // analysisThreads()'s choice; results do not depend on it).
    std::size_t threads() const { return m_threads; }
    void setThreads(std::size_t threads) { m_threads = threads; }
    // instance -> cache key of the last analyze() (dict order)
    const std::vector<std::pair<std::string, std::string>>& analysisKeys() const { return m_keys; }

    // check(): PipelineError("Invalid parameters:\n  - a\n  - b") when core::validate finds
    // problems, then "WARNING: <text>" for every layout warning.
    void check() const;

    // resolve_inputs(reporter, force_clean): the input files, with cleaning the cleaned ones
    // ("Cleaning trajectory with cpptraj ...", then analysis::clean with this reporter; `forceClean`
    // runs cpptraj again), then progress(1, "inputs ready"). Kept for later calls unless
    // forceClean. Throws CpptrajError (cleaning.hpp), core::Cancelled.
    const ResolvedInputs& resolveInputs(const core::Reporter* reporter = nullptr, bool forceClean = false);
    // cleaned_if_available: the cleaned files when a finished cleaning with these inputs and
    // settings exists (its marker naming this key, MDATK 0.1.0), else the raw inputs; never runs
    // cpptraj (previews, the interface's structure view).
    ResolvedInputs cleanedIfAvailable() const;
    // universe(): the loaded system of the resolved inputs (kept). Throws what io::loadSystem
    // throws, and core::ValueError("No trajectory loaded into Universe") without coordinates. The
    // readers' warnings (core::System::warnings) go to the reporter as "WARNING: <message>" once,
    // when the system is loaded.
    core::System& system();

    // input_fingerprint: {"topology": fileFingerprint, "trajectories": [...], "frames": {"start",
    // "stop", "step"}, "dt_ps"}, plus "time_axis": "frame x dt_ps" when inputs.dt_ps overrides the
    // time step (MDATK 0.1.0; what every cache key covers).
    core::Value inputFingerprint(const ResolvedInputs& inputs) const;

    // required_analyses(video, charts): the enabled user analyses, the panels' analyses when the
    // video is on, the instances of charted series ("<instance>/<output>") when charts are on
    // (video / charts: the project's output settings when not given), without "external";
    // instances instanceProject rejects are left out with the log line
    // "WARNING: skipping data source '<instance>': <error>" (a KeyError's message without the
    // quotes MDATK prints, CLI-01; in sorted order; MDATK's order is
    // that of a Python set). Sorted by registry index, then instance id.
    std::vector<std::string> requiredAnalyses(std::optional<bool> video = std::nullopt,
                                              std::optional<bool> charts = std::nullopt) const;

    // analyze(names, force, reporter): `names` (requiredAnalyses() when not given) are loaded from
    // the cache ("<label>: loaded from cache (<instance>-<key>)") unless `force`, the missing ones
    // run in one pass over the frames of inputs.start / stop / step (runAnalyses; the reporter's
    // 0.02 - 1 part) and are stored; then the trajectory is back on frame 0, as MDATK's iteration
    // over a frame list leaves it. Without names, the frame axis alone (times read frame by frame).
    // Then the dataset (time unit inputs.time_unit), display overrides, external data
    // (PipelineError("external data: <reason>"); the files are first read once before the frame
    // loop, so such an error comes before the analyses run: CLI-02), the manifest, progress(1, "analysis done").
    // Throws std::out_of_range / core::ValueError for unknown analyses and bad instances,
    // core::ValueError for empty frame selections, core::Cancelled, and what the analyses and the
    // readers throw.
    core::Dataset analyze(const std::optional<std::vector<std::string>>& names = std::nullopt, bool force = false,
                          const core::Reporter* reporter = nullptr);

    // apply_display_overrides: analyses.energy.target_value as the "target" of every energy
    // instance's "total" series (null when not set).
    void applyDisplayOverrides(core::Dataset& dataset) const;

    // load_cached_dataset: the dataset of the latest run, from the manifest (display overrides
    // applied). Throws what CacheStore::loadDatasetFromManifest throws.
    core::Dataset loadCachedDataset() const;

    // CIT-06 (docs/DIFFERENCES.md): the citation keys the cached results of the last analyze()
    // record (meta "citation_keys": what a run used, e.g. the energy analysis' force fields, the
    // covalent radii of a bond guess), in analysisKeys() order, each result's keys in its order,
    // without repeats. Entries missing from the cache, or without the key, add nothing; empty
    // before analyze().
    std::vector<std::string> runCitationKeys() const;

private:
    core::Project m_project;
    core::Reporter m_reporter;
    CacheStore m_cache;
    std::optional<ResolvedInputs> m_inputs;
    std::unique_ptr<core::System> m_system;
    std::vector<std::pair<std::string, std::string>> m_keys;
    std::size_t m_threads = 0;
};

}  // namespace mdbridge::analysis
