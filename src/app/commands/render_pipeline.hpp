#pragma once
// The output part of MDATK's end-to-end pipeline (mdatk/pipeline.py Pipeline: check, stamp,
// render_video, render_charts, export_pdbs, render_preview, write_params_log, run) on MDBridge's
// analysis pipeline (mdbridge/analysis/pipeline.hpp), for the `render`, `charts` and `preview`
// commands (and the GUI later). Qt-free: what draws and encodes is a RenderBackend, which the
// `mdbridge` program fills with its Qt implementation (app/render_backend_qt.cpp) and the tests with
// their own.
//
// Kept from MDATK: the run's order and progress ranges (inputs 0-0.15, analyses 0.15-0.55, video
// 0.55-0.92, charts 0.92-0.97, PDBs 0.97-1), one time stamp for every file of a run, the file names
// (output.video.filename stamped, charts, PDBs, MDATK_params_log.toml), the log lines ("<program>
// <version> - output directory D", "Finished in N s"), the params log's run table and the step log
// entry "render" / "Analyses and outputs" (files, summary "N frames analysed; M file(s) written to
// D", the pipeline's citations). Differences (docs/DIFFERENCES.md "Rendering: the render, charts and
// preview commands"):
//
//   RN-02  no analysed frames: "no analysed frames to render" (PipelineError) before any figure is
//          built, for the video, the charts and the preview
//   RN-11 ... RN-15  check() adds render::panels::renderProblems to core::validate's problems; the
//          runs from the cache (which skip check(), as MDATK) check those alone
//   D8     output.video.ffmpeg_path is not used (the bundled FFmpeg libraries encode); a project that
//          names another executable gets the log line "Note: output.video.ffmpeg_path (P) is not
//          used; MDBridge encodes with its bundled FFmpeg libraries"
//   (milestone 6: a shown density map without a placement is fitted first, MDATK's ensure_density_fit,
//          density_pipeline.hpp; milestone 4's D9 warning is gone)
//   environment  the params log's environment table is MDBridge's (program, platform, the
//          backend's libraries), not Python's packages and `ffmpeg -version`

#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/text.hpp"

namespace mdbridge::app {

// What draws, writes and encodes figures. Every member may be empty: drawing then fails with
// "figures cannot be drawn by this program" (a build without Qt), the video with "video encoding
// is not available in this build".
struct RenderBackend {
    // Called once before anything is drawn (the Qt application on the offscreen platform); the
    // object is kept until the command ends.
    std::function<std::shared_ptr<void>()> start;
    // The figure fonts (the bundled IBM Plex Serif faces and their fallback).
    std::function<const render::FontSet&()> fonts;
    // Writes a drawn figure (render::Figure::draw) to `path` in `format` (one of
    // render::charts::chartFormats()); throws naming the file when it cannot.
    std::function<void(const render::Figure& figure, const std::filesystem::path& path, const std::string& format)>
        writeFigure;
    // MDATK's render_video(dataset, project, path, reporter): returns the written file; throws
    // (VideoError) on failure, core::Cancelled when cancelled. The video frame loop of the program.
    std::function<std::filesystem::path(const core::Dataset& dataset, const core::Project& project,
                                        const std::filesystem::path& path, const core::Reporter& reporter)>
        renderVideo;
    // Entries added to the params log's environment table (library versions).
    std::function<core::Table()> environment;
};

// RunOutputs: what a run wrote.
struct RunOutputs {
    std::optional<std::string> video;
    std::vector<std::string> charts;
    std::vector<std::string> pdbs;
    std::optional<std::string> manifest;
    std::optional<std::string> paramsLog;
    std::optional<core::Dataset> dataset;

    // summary(): {"video", "charts", "pdbs", "cache_manifest", "params_log"} (null when not written)
    core::Value summary() const;
};

// stamped("MDATK_video.mp4", stamp) = "MDATK_video_<stamp>.mp4" ("" stamp: unchanged).
std::string stamped(const std::string& name, const std::string& stamp);

// check(): PipelineError("Invalid parameters:\n  - a\n  - b") for the problems of core::validate
// and render::panels::renderProblems, then "WARNING: <text>" for every layout warning.
void checkProject(const core::Project& project, const core::Reporter& reporter);
// The render problems alone (runs from the cache): PipelineError as above.
void checkRenderSettings(const core::Project& project);

class RenderPipeline {
public:
    // `projectFile`: the project file whose step log records a run (empty: none, as for a params
    // log); `backend` must outlive the pipeline.
    RenderPipeline(analysis::Pipeline& analyses, const RenderBackend& backend, std::filesystem::path projectFile = {});

    analysis::Pipeline& analyses() { return m_analyses; }
    const core::Project& project() const { return m_analyses.project(); }

    // check() (checkProject with the analysis pipeline's reporter).
    void check() const;

    // stamp(): the run's time stamp ("" when output.timestamp_files is off), made once; newStamp()
    // forgets it.
    std::string stamp();
    void newStamp() { m_stamp.reset(); }
    // output_path(name): `name` when absolute, else out_dir / name.
    std::filesystem::path outputPath(const std::string& name) const;

    // render_video: RN-02, D8 note, ensure_density_fit, the backend's renderVideo on
    // out_dir / stamped(output.video.filename). Returns the path written.
    std::string renderVideo(const core::Dataset& dataset, const core::Reporter* reporter = nullptr);
    // render_charts: RN-02, render::charts::renderCharts into out_dir with the run's stamp.
    std::vector<std::string> renderCharts(const core::Dataset& dataset, const core::Reporter* reporter = nullptr);
    // export_pdbs: output.pdb's frames / times / selection of the analysis pipeline's system.
    std::vector<std::string> exportPdbs(const core::Reporter* reporter = nullptr);
    // render_preview: ensure_density_fit, RN-02, the figure at output.video's size and dpi, panels updated
    // to `pos`, written to `path` (format from its extension; none: png, written as path + ".png").
    // Returns the path written.
    std::string renderPreview(const core::Dataset& dataset, std::size_t pos, const std::string& path);
    // write_params_log: <out_dir>/MDATK_params_log.toml. Returns its path.
    std::string writeParamsLog(const RunOutputs& outputs);

    // run(video, charts, pdb, force_analysis): check(), then the run inside the step log entry
    // "render" (std::nullopt: the project's output.*.enabled).
    RunOutputs run(std::optional<bool> video, std::optional<bool> charts, std::optional<bool> pdb, bool forceAnalysis);

private:
    void startBackend();
    const render::FontSet& fonts();
    RunOutputs runSteps(bool video, bool charts, bool pdb, bool forceAnalysis);

    analysis::Pipeline& m_analyses;
    const RenderBackend& m_backend;
    std::filesystem::path m_projectFile;
    std::optional<std::string> m_stamp;
    std::shared_ptr<void> m_started;
    bool m_startCalled = false;
    bool m_inputsResolved = false;
};

// The program's own entries of the params log's environment table: {"mdbridge": version,
// "platform": "<system> <machine>"} (the backend adds its libraries).
core::Table environmentInfo();

}  // namespace mdbridge::app
