#pragma once
// The analysis window's jobs and desktop actions (milestone 5 W2-B; MDATK v0.0.4.1 gui/app.py
// render, run_cleaning, inspect_system, preview_frame, export_pdbs, rerender_from_cache,
// open_output_dir, open_video, centre_density, fit_density, start_job; gui/desktop.py latest_video).
//
// Every job runs on the window's jobs::JobRunner (one at a time, progress, log, Cancel, the close
// question AP-14) through the analysis and render pipelines the command line uses
// (analysis::Pipeline, RenderPipeline with the Qt RenderBackend). Before a job the typed text is
// stored (FormContext::commitPendingEdits) and, as v0.0.4.1's collect() (owner decision D7), text
// a field could not take blocks every job; Render, Render still and Export PDBs are also blocked by
// the problems of core::validate and render::panels::renderProblems (what the command line's
// check() refuses), with MDATK's titles: "Please fix these parameters before rendering", "Please
// fix these parameters first", "Invalid parameters". The result dialogs are MDATK's ("Render
// finished", "Cleaning finished" / "Cleaned files reused", "PDB export finished", "Re-render
// finished"); Render still shows the frame in an ImageDialog.
//
// Differences from MDATK (docs/DIFFERENCES.md WA-01 ...): the job buttons are disabled while a
// job runs (MDATK: Run cleaning now and Inspect system answered "Another job is still running.");
// Open video prefers the video of the last render or re-render only while it is in the current
// output folder (MDATK: also after the folder was changed, and a re-render did not count); Inspect
// system unfolds the log, where its result is, and reports an empty group as such; after Run
// cleaning the reference frame is loaded again from the cleaned files; the folder and the video
// open through QDesktopServices.
//
// Global search and Refine (milestone 6 W3-A) run the density fit of the command line
// (density_pipeline.hpp fitDensity: chem::fitModel, progress and Cancel through the job's
// core::Reporter) as the job "Fitting the density map", blocked as v0.0.4.1's fit_density is
// (collect(): "Please fix these parameters first"; no map: "No density map"). The placement found
// (transform, anchor, fit_cc, fit_info: MDATK's _store_density_result) is stored as one undo step
// once the job has succeeded, so a cancelled or failed fit leaves the project as it was (no undo
// step, no modified marker); a map file changed during the fit gets no placement of the old one
// (WF-02). The preview is drawn again with the placement, and the reference frame loaded again
// when the fit cleaned the inputs anew (MDATK: pipe.structure_preview(); WF-03). Render and Render
// still keep the placement an automatic fit found (density.auto_fit, MDATK's _adopt_auto_fit:
// "Density map placement stored: <fit_info>"), also as one undo step.
//
// The pages' buttons are connected here: Input & cleaning (runCleaningRequested,
// inspectSystemRequested), Outputs (exportPdbsRequested, rerenderRequested), Structure view
// (fitRequested, centreOnModelRequested) and the Layout page's *Render still*: a page with a signal
// `renderStillRequested(int frame)` is connected by name, and its button named
// "layout.render_still" (if any) is disabled while a job runs.

#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>

#include <QObject>
#include <QPointer>
#include <QString>
#include <QUrl>

#include "mdbridge/core/project.hpp"

namespace mdbridge::app {

class ImageDialog;
class MainWindow;
struct RenderBackend;

class WindowActions : public QObject {
    Q_OBJECT

public:
    explicit WindowActions(MainWindow& window);
    ~WindowActions() override;
    WindowActions(const WindowActions&) = delete;
    WindowActions& operator=(const WindowActions&) = delete;

    // ---- the actions (each returns false when nothing was started: refused, or a job runs) ----

    bool render();                 // Render (Ctrl+R): cleaning, analyses, every enabled output
    bool renderStill(int frame);   // Render still: one frame of the video, shown in a dialog
    bool exportPdbs();             // Export PDBs now
    bool rerenderFromCache();      // Re-render video + charts from cache
    bool runCleaning();            // Run cleaning now
    bool inspectSystem();          // Inspect system
    void openOutputFolder();       // Open output folder
    void openVideo();              // Open video
    void centreOnModel();          // Structure view › Density map › Centre on model
    // Structure view › Density map › Global search ("global") / Refine ("manual")
    bool fitDensity(const QString& search);

    // ---- state ----

    // The video the last Render or Re-render wrote ("" before one).
    const std::filesystem::path& lastVideo() const { return m_lastVideo; }
    // The video Open video opens now ("" when there is none).
    std::filesystem::path videoToOpen();
    // The still dialog shown last (null once closed).
    ImageDialog* stillDialog() const;
    // Centre on model can work: the reference frame is loaded and a map file is chosen.
    bool centreAvailable() const;

    // How folders and files are opened (QDesktopServices::openUrl by default; tests record).
    using UrlOpener = std::function<bool(const QUrl& url)>;
    void setUrlOpener(UrlOpener opener);

    // MDATK's latest_video(out_dir, filename): the newest file in `outDir` whose name is the
    // base name of `filename` with anything between its stem and extension (the time stamp), not
    // hidden (the encoder's temporary .MDATK_tmp_* files); equal times: the last name. "" when none.
    static std::filesystem::path latestVideo(const std::filesystem::path& outDir, const std::string& filename);
    // The still's file: <temporary folder>/MDATK_still_<process id>.png (MDATK's name).
    static std::filesystem::path stillPath();

signals:
    // A job of these actions finished with its result dialog shown (`title`: "Render finished" ...).
    void resultShown(const QString& title);
    // The still dialog was opened.
    void stillShown(const QString& path);

private slots:
    void onRenderStill(int frame);  // the Layout page's renderStillRequested(int)

private:
    void connectPages();
    void updateCentre();
    // The placement of a fit (transform, anchor, fit_cc, fit_info) stored as the undo step `label`
    // when it is for the map the settings name now; false (and a warning in the log) otherwise.
    bool storeFit(const core::DensityParams& fitted, const std::string& label);
    // MDATK's _adopt_auto_fit: a placement an automatic fit of a render job found.
    void adoptAutoFit(const core::DensityParams& fitted);
    // MDATK's collect(validate): the project (paths resolved) or nothing after the error dialog.
    std::optional<core::Project> collect(bool validate, bool renderChecks, const QString& title);
    void result(const QString& title, const QString& detail);
    std::filesystem::path projectFile() const;

    MainWindow& m_window;
    std::shared_ptr<RenderBackend> m_backend;
    UrlOpener m_opener;
    std::filesystem::path m_lastVideo;
    QPointer<ImageDialog> m_still;
};

}  // namespace mdbridge::app
