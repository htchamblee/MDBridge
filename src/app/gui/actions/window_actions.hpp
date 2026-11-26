#pragma once
// The analysis window's jobs and desktop actions (milestone 5 W2-B; MDATK v0.0.4.1 gui/app.py
// render, run_cleaning, inspect_system, preview_frame, export_pdbs, rerender_from_cache,
// open_output_dir, open_video, centre_density, start_job; gui/desktop.py latest_video).
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
// The pages' buttons are connected here: Input & cleaning (runCleaningRequested,
// inspectSystemRequested), Outputs (exportPdbsRequested, rerenderRequested), Structure view
// (centreOnModelRequested) and the Layout page's *Render still*: a page with a signal
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
