// `mdbridge` (console: the command line and `mdbridge gui`) and, on Windows, `mdbridge-gui`
// (no console window: the interface).
#include "cli.hpp"
#include "mdbridge/render/panels/vulkan_stage.hpp"

int main(int argc, char** argv) {
    // The structure view's 3D layer (figures, video frames, previews) on the GPU when Vulkan can run
    // here, else the built-in software renderer (render/panels/vulkan_stage.hpp). The device opens
    // on first use only, and is released when main returns, before static destruction.
    const auto gpu = mdbridge::render::panels::installGpuStageRenderers();
#ifdef MDBRIDGE_GUI_ENTRY
    return mdbridge::app::guiEntryMain(argc, argv);
#else
    return mdbridge::app::cliMain(argc, argv);
#endif
}
