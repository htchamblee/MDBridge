#pragma once
// `mdbridge gui PROJECT --script FILE`: a short list of user actions run against the window, for
// the smoke test of finished builds (packaging/smoke_test.py --gui) and for support. The actions
// go through the window as a user's would: values are typed into the field widgets (key events,
// Enter), check boxes toggled with Space, Save and Render still are the window's own actions.
//
// The file is JSON:
//
//   {"format": "mdbridge-gui-script", "version": 1, "actions": [
//     {"action": "wait"},                                     project loaded, preview drawn
//     {"action": "page", "page": "outputs"},                  show a page (core::kProblemPages ids)
//     {"action": "set", "key": "output.video.fps", "value": 30},   enter a value (the field's page is shown)
//     {"action": "type", "key": "notes", "text": "..."},      type text into a field and commit it
//     {"action": "save"},                                      Save (Ctrl+S)
//     {"action": "copy_project", "to": "before.toml"},        copy the project file as it is on disk
//     {"action": "render_still", "frame": 0, "to": "still.png", "screenshot": "still"},
//                                                              Render still: its image copied, the
//                                                              dialog photographed, then closed
//     {"action": "screenshot", "name": "window"},             the window as <output>/window.png
//     {"action": "screenshot", "name": "page", "pages": true} every page: <output>/page-<id>.png
//   ]}
//
// Relative paths are taken from the script's folder; screenshots go to `output` (default: the
// script's folder). Questions are answered by the action's "answers" ({"kind": "choice id"}), else
// with the question's safe default; an error message fails the action. A JSON report is printed on
// one line at the end: {"ok", "problem", "project", "modified", "pages", "actions": [{"action", "ok",
// ...}], "messages", "questions", "libraries" (the shared libraries loaded, Linux: which copies an
// AppImage used)}; the window then quits with 0 when every action succeeded, else 1 (2: the script
// could not be read). The report's paths ("project", "file", "files", "image", "dialog") are
// stored paths (as the project file and settings keep them, '/' on Windows); its texts ("problem",
// "messages", "questions") and stderr show paths as the interface does (gui::displayPath: '\' on
// Windows, docs/DIFFERENCES.md WN-12).

#include <filesystem>
#include <memory>
#include <string>

namespace mdbridge::app {

class MainWindow;

class ScriptRunner {
public:
    // Reads the script (throws std::runtime_error with the reason when it cannot be read).
    ScriptRunner(MainWindow& window, const std::filesystem::path& script, const std::filesystem::path& output);
    ~ScriptRunner();
    ScriptRunner(const ScriptRunner&) = delete;
    ScriptRunner& operator=(const ScriptRunner&) = delete;

    // Runs every action (waiting inside, as the GUI tests do), prints the report and returns the
    // exit code: 0 when every action succeeded, 1 when one failed (the rest are not run).
    int run();

    // After run(): the report's "problem" ("" when every action succeeded).
    const std::string& problem() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

}  // namespace mdbridge::app
