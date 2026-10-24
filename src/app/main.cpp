// `mdbridge` (console: the command line and `mdbridge gui`) and, on Windows, `mdbridge-gui`
// (no console window: the interface).
#include "cli.hpp"

int main(int argc, char** argv) {
#ifdef MDBRIDGE_GUI_ENTRY
    return mdbridge::app::guiEntryMain(argc, argv);
#else
    return mdbridge::app::cliMain(argc, argv);
#endif
}
