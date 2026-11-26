#pragma once
// What the commands share (commands.cpp and the chem commands of chem_commands.cpp,
// density_commands.cpp): MDATK's _common() options, load_project(), the console reporter
// (make_reporter), the step log's command text. Internal to mdbridge::commands.

#include <chrono>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "commands.hpp"
#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::app::detail {

// The options MDATK's _common() adds to init, info, export-pdb, ... density fit / simulate.
struct CommonArguments {
    std::string params;
    std::string top;
    std::vector<std::string> traj;
    std::string frames;
    bool clean = false;
    bool noClean = false;
    std::string cpptraj;
    std::string condaEnv;
    std::string outDir;
    std::string cacheDir;
    std::vector<std::string> set;
    bool quiet = false;
    bool verbose = false;
    CLI::Option* cleanOption = nullptr;
    CLI::Option* noCleanOption = nullptr;
};


// A single-value option given more than once takes the last value, as argparse does.
CLI::Option* lastWins(CLI::Option* option);
// _common(p)
void addCommon(CLI::App& command, CommonArguments& a);
// os.path.abspath
std::string absolutePath(const std::string& path);
// load_project(args, resolve)
core::Project loadProject(const CLI::App& command, const CommonArguments& a, bool resolve);

// The console reporter of MDATK's command line (make_reporter, ConsoleProgress).
class ConsoleReporter {
public:
    ConsoleReporter(Console& console, bool quiet, bool verbose)
        : m_console(console), m_quiet(quiet), m_verbose(verbose) {}

    core::Reporter reporter() {
        return core::Reporter([this](const std::string& text) { logLine(text); },
                              [this](double fraction, const std::string& message) { progress(fraction, message); });
    }

private:
    void logLine(const std::string& text) {
        const bool detail = text.starts_with("    ") || text.starts_with("\t") || text.starts_with("|");
        if ((m_verbose || !detail) && !m_quiet) {
            m_console.err << (m_console.errIsTerminal ? "\r\033[K" : "") << text << '\n';
            m_console.err.flush();
        }
    }

    void progress(double fraction, const std::string& message) {
        if (m_quiet) {
            return;
        }
        const double now =
            std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
        if (fraction < 1.0 && fraction - m_last < 0.01 && now - m_lastTime < 2.0) {
            return;
        }
        m_last = fraction;
        m_lastTime = now;
        char percent[32];
        std::snprintf(percent, sizeof percent, "%5.1f", fraction * 100);
        const std::string line = "[" + std::string(percent) + "%] " + message;
        if (m_console.errIsTerminal) {
            m_console.err << "\r\033[K" << line << (fraction >= 1.0 ? "\n" : "");
        } else if (fraction >= 1.0) {
            // MDATK also prints when int(fraction * 10) != int(self.last * 10), compared after
            // self.last = fraction: never, so only final lines reach a file or a pipe
            m_console.err << line << '\n';
        }
        m_console.err.flush();
    }

    Console& m_console;
    bool m_quiet;
    bool m_verbose;
    double m_last = -1.0;
    double m_lastTime = 0.0;
};

// The step log's command text: " ".join(sys.argv[1:]) when the command line's first word is
// `group` (MDATK: entry.command = ... if sys.argv[1:2] == [group] else "").
std::string commandText(const std::vector<std::string>& arguments, const std::string& group);

}  // namespace mdbridge::app::detail
