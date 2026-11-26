#pragma once
// The commands of MDATK's `mdatk` command line that the libraries support so far: demo, init
// (init-params), info, list, clean, analyze, export-pdb, render, charts, preview, density info /
// fit / simulate, the chemistry commands pdb, build and dict (milestone 6: chem_commands.hpp), and
// MDBridge's citations (D14: the references of a project's settings and cached results).
// Names, options, defaults, output, files and exit codes are MDATK v0.0.4.1's (mdatk/cli.py);
// MDATK's workflow commands follow with milestone 7 and are not offered until then. render, charts, preview and demo --run draw with a
// RenderBackend (render_pipeline.hpp) that the program passes (the Qt one in the `mdbridge`
// program); without one they fail with "figures cannot be drawn by this program". `analyze` keeps
// MDBridge's own cache (mdbridge/analysis/cache.hpp); cleaning runs cpptraj as MDATK does
// (mdbridge/analysis/cleaning.hpp).
//
// This library needs no Qt: the `mdbridge` program adds these commands to its own (gui, project,
// graphics), and the tests run them in-process.

#include <exception>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/project.hpp"
#include "pubchem.hpp"

namespace mdbridge::app {
struct RenderBackend;
namespace detail {
class ChemCommands;
}
}

namespace CLI {
class App;
}

namespace mdbridge::app {

// Where a command writes: its result to `out`, messages and progress to `err` (MDATK's stdout and
// stderr). `errIsTerminal`: progress lines are redrawn in place (\r, ESC[K) as on a terminal.
struct Console {
    std::ostream& out;
    std::ostream& err;
    bool errIsTerminal = false;
};

// The process's standard streams (errIsTerminal from isatty).
Console standardConsole();

class CoreCommands {
public:
    // Adds the commands to `app` (the root of the command line). `network`: how `dict pubchem`
    // reaches PubChem (the `mdbridge` program's Qt Network; null: no lookups). Both must outlive
    // the commands.
    explicit CoreCommands(CLI::App& app, const RenderBackend* backend = nullptr,
                          const NetworkBackend* network = nullptr);
    ~CoreCommands();
    CoreCommands(const CoreCommands&) = delete;
    CoreCommands& operator=(const CoreCommands&) = delete;

    // After app.parse(): runs the command that was given and returns its exit status (0, 1 for
    // errors, which are printed as "error: <message>", 2 for usage problems, 130 when cancelled);
    // std::nullopt when none of these commands was given. Warnings logged meanwhile (core::logWarning)
    // go to console.err as "WARNING: <message>".
    std::optional<int> run(Console& console);

    // The command line's words after the program name (UTF-8), for the step log's `command`
    // (MDATK: " ".join(sys.argv[1:])).
    void setArguments(std::vector<std::string> arguments) { m_arguments = std::move(arguments); }

    struct Arguments;  // the parsed options (commands.cpp)

private:
    std::unique_ptr<Arguments> m_args;
    std::unique_ptr<detail::ChemCommands> m_chem;
    std::vector<std::string> m_arguments;
};

// A command line with only these commands (the program name is argv[0]), for tests and for
// checking parity without the Qt front end: parse errors print CLI11's message and give 2, as
// argparse does; without a command the help is printed and 2 returned (MDATK's `mdatk`).
int runCoreCommandLine(int argc, const char* const* argv, Console& console, const RenderBackend* backend = nullptr,
                       const NetworkBackend* network = nullptr);

// What `list` prints (cmd_list): the analyses, the conda environments with cpptraj
// (core::conda::findCpptrajEnvs), the panel types, the layout presets and a note on data keys.
void printList(std::ostream& out);

// MDATK's apply_override(project, "dotted.path=value"): the value is read as JSON when it is
// valid JSON (else it is the text), list elements are addressed by index (groups.0.residues),
// a true/false setting needs a JSON boolean, a text setting takes the text as given; the result is
// read through projectFromPlain (its type checks and warnings). Throws core::ValueError with
// MDATK's messages ("--set expects key=value, got 'x'", "unknown parameter 'a.b'", "a.b expects
// true/false", "groups.5: no list element '5'", "empty key in '=1"), core::ProjectError.
void applyOverride(core::Project& project, std::string_view assignment);

// Python's json.dumps(value, indent=2) (ensure_ascii: non-ASCII characters and DEL as \uXXXX
// escapes, surrogate pairs beyond U+FFFF, bytes that are not UTF-8 as \udcXX): what MDATK prints.
std::string jsonForConsole(const core::Value& value);

// str(exc) of the Python exception MDATK raises where C++ throws `error` (OSError texts for
// file-system errors, the message otherwise). A KeyError (std::out_of_range) gives its message
// without the quotes of Python's str(KeyError) (docs/DIFFERENCES.md, CLI-01).
std::string pythonErrorText(const std::exception& error);

}  // namespace mdbridge::app
