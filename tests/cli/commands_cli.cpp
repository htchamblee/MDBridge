// `mdbridge_commands_cli`: the MDATK commands of the command line (init, info, export-pdb,
// density info) without the Qt front end, for tests/parity/cli_parity.py on builds without Qt
// (MDBRIDGE_BUILD_APP=OFF). The `mdbridge` program offers the same commands.
#include "commands.hpp"

int main(int argc, char** argv) {
    mdbridge::app::Console console = mdbridge::app::standardConsole();
    return mdbridge::app::runCoreCommandLine(argc, argv, console);
}
