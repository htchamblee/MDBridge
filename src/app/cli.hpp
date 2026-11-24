#pragma once
// The command line front end: `mdbridge COMMAND ...`. Command names and options follow MDATK's
// `mdatk` command as they are ported: gui, project and graphics (milestone 1); init, info,
// export-pdb and density info (milestone 2, commands/commands.hpp); the rest with later milestones.

namespace mdbridge::app {

// `mdbridge` (console program).
int cliMain(int argc, char** argv);

// `mdbridge-gui` (Windows program without a console): the options of `mdbridge gui`.
int guiEntryMain(int argc, char** argv);

}  // namespace mdbridge::app
