#pragma once
// MDATK's chemistry commands on the chem library (milestone 6, docs/specs/m6-plan.md 1.6):
//   pdb clean / check / to-lrna / annotate    (mdatk/pdbtools.py run, cli.py cmd_pdb, _pdb_annotate)
//   build helix / mutate / fit                (mdatk/builder.py run, cli.py cmd_build)
//   dict list / smiles / finish / pubchem     (mdatk/dictionaries.py run, cli.py cmd_dict)
//   density fit / simulate                    (cli.py cmd_density_fit, cmd_density_simulate)
// with MDATK's options, defaults, output (notes on stderr prefixed "pdbtools:", "annotate:",
// "build:", "dictionaries:"), files and the step-log entries of `--project` (step, title, command,
// files, summary, citations). MDBridge's differences are those of the library (CH-xx, D7, D13) and
// of the command line: CH-10 (resolutions and spacings), CH-12 (`density fit --save` through the
// edit model), CH-13 (a missing input file is named). Internal to mdbridge::commands.

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "commands.hpp"
#include "pubchem.hpp"

namespace CLI {
class App;
}

namespace mdbridge::app::detail {

class ChemCommands {
public:
    // Adds pdb, build and dict to `app` and fit / simulate to `density`. `network` may be null (no
    // PubChem lookups).
    ChemCommands(CLI::App& app, CLI::App& density, const NetworkBackend* network);
    ~ChemCommands();
    ChemCommands(const ChemCommands&) = delete;
    ChemCommands& operator=(const ChemCommands&) = delete;

    // One of these commands (or its group alone) was given.
    bool parsed() const;
    // density fit / simulate was given (density info is the core commands').
    bool densityParsed() const;
    // Runs the command given; `arguments` are the command line's words after the program name (the
    // step log's command). Errors propagate (CoreCommands::run reports them); a group without its
    // action prints the group's help and gives 2.
    int run(Console& console, const std::vector<std::string>& arguments);

    struct Arguments;

private:
    std::unique_ptr<Arguments> m_args;
};

}  // namespace mdbridge::app::detail
