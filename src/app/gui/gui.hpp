#pragma once
// The Qt interface.

#include <string>

#include "mdbridge/render/view.hpp"

namespace mdbridge::app {

struct GuiOptions {
    std::string project;  // project file to open (UTF-8), "" for none
    render::Renderer renderer = render::Renderer::Auto;
    // Checks for CI and support: save the 3D view after `screenshotFrames` frames to this PNG,
    // print a JSON report and quit.
    std::string screenshot;
    int screenshotFrames = 10;
    double quitAfter = 0;  // quit after this many seconds (0: run until closed)
};

// Runs the interface; `argc` / `argv` go to QApplication (Qt's own options such as -platform).
int runGui(const GuiOptions& options, int& argc, char** argv);

}  // namespace mdbridge::app
