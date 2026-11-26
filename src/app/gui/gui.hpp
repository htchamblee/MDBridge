#pragma once
// The Qt interface.

#include <string>

#include "mdbridge/render/view.hpp"

namespace mdbridge::app {

struct GuiOptions {
    std::string project;  // project file to open (UTF-8), "" for none
    render::Renderer renderer = render::Renderer::Auto;
    // Checks for CI and support: save the 3D view to this PNG once the project's structure and
    // figure preview are loaded (and the view shows the structure) and at least `screenshotFrames`
    // frames were drawn, print a JSON report and quit; `previewScreenshot`: the preview too.
    std::string screenshot;
    std::string previewScreenshot;
    int screenshotFrames = 10;
    double quitAfter = 0;  // quit after this many seconds (0: run until closed)
    // Check: run this action script (gui/script_runner.hpp) once the window is shown, print its
    // report and quit (0: every action succeeded); screenshots into `scriptOutput` ("" = the
    // script's folder).
    std::string script;
    std::string scriptOutput;
};

// Runs the interface; `argc` / `argv` go to QApplication (Qt's own options such as -platform).
int runGui(const GuiOptions& options, int& argc, char** argv);

}  // namespace mdbridge::app
