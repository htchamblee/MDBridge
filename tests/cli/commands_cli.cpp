// `mdbridge_commands_cli`: the MDATK commands of the command line (app::runCoreCommandLine, every
// command of commands.cpp: demo, init, info, list, clean, analyze, export-pdb, render, charts,
// preview, citations, density, pdb, build, dict) without the Qt front end, for
// tests/parity/cli_parity.py on builds without Qt (MDBRIDGE_BUILD_APP=OFF); the commands that draw or
// fetch report what the program lacks. The `mdbridge` program offers the same commands.
#include "commands.hpp"

int main(int argc, char** argv) {
    mdbridge::app::Console console = mdbridge::app::standardConsole();
    return mdbridge::app::runCoreCommandLine(argc, argv, console);
}
