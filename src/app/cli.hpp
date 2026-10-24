#pragma once
// The command line front end: `mdbridge COMMAND ...`. Command names and options follow MDATK's
// `mdatk` command as they are ported (milestone 5); milestone 1 has gui, project and graphics.

namespace mdbridge::app {

// `mdbridge` (console program).
int cliMain(int argc, char** argv);

// `mdbridge-gui` (Windows program without a console): the options of `mdbridge gui`.
int guiEntryMain(int argc, char** argv);

}  // namespace mdbridge::app
