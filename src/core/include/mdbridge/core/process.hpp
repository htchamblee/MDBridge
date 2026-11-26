#pragma once
// Running external programs (cpptraj for the trajectory cleaning; the planner's scripts in
// milestone 7): finding a program on PATH as shutil.which does, and a process runner with merged
// stdout / stderr delivered as it arrives, cancellation that does not wait for the program to
// print, and the program's whole process tree ended on cancel. Qt-free; Linux (fork / exec,
// process groups) and Windows (CreateProcessW, a job object).
//
// MDATK runs cpptraj with subprocess.Popen and a blocking read loop (mdatk/cleaning/cpptraj.py
// run_cpptraj): a cancel there is noticed only when cpptraj next prints, and `conda run`'s child
// survives terminate(). See docs/DIFFERENCES.md (CLN-03).

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mdbridge::core::process {

// A program that could not be started (not found, not executable, no such working directory):
// what() is Python's OSError text ("[Errno 2] No such file or directory: 'x'").
class StartError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// os.environ.get(name): "" when unset (UTF-8 on Windows).
std::string environmentVariable(std::string_view name);

// os.access(path, X_OK) and a regular file (Windows: an existing file).
bool isExecutableFile(const std::string& path);

// shutil.which(name): a name with a directory part is checked as it is; otherwise the folders of
// PATH (os.defpath "/bin:/usr/bin" without PATH) in order, the first executable regular file
// (Windows: the current folder first, and the extensions of PATHEXT when the name has none of
// them). std::nullopt when none is found.
std::optional<std::string> which(std::string_view name);

// The environment for programs MDBridge starts: the process's own, without what the AppImage's
// start script added (LD_LIBRARY_PATH / PATH entries inside $APPDIR, FONTCONFIG_FILE inside it),
// as MDATK's frozen.child_env drops its bundle's variables. Pairs NAME=value in the process's order.
std::vector<std::pair<std::string, std::string>> childEnvironment();

struct RunOptions {
    // Working directory ("" = the current one).
    std::filesystem::path workingDirectory;
    // The environment (std::nullopt: childEnvironment()).
    std::optional<std::vector<std::pair<std::string, std::string>>> environment;
    // Output as it arrives (stdout and stderr merged, bytes as the program wrote them).
    std::function<void(std::string_view chunk)> onOutput;
    // Polled at least every `pollInterval`; true ends the program (see run()).
    std::function<bool()> cancelled;
    std::chrono::milliseconds pollInterval{50};
    // After a cancel: time between asking the program to stop (SIGTERM to its process group) and
    // killing it (SIGKILL). Windows ends the job's processes at once (TerminateJobObject), as
    // Python's terminate() does there.
    std::chrono::milliseconds terminateGrace{10000};
};

struct Result {
    // Exit status; a program ended by a signal gives -signal (Python's returncode).
    int exitCode = 0;
    bool cancelled = false;  // ended because options.cancelled() returned true
};

// Runs argv[0] (a path, or a name looked up with which()) with the arguments argv[1..], stdin
// from the null device. Returns when the program and the output pipe are done, or after a cancel
// once the program has ended. Throws StartError when it cannot start, std::invalid_argument for an
// empty argv.
Result run(const std::vector<std::string>& argv, const RunOptions& options = {});

// Splits output chunks into lines as MDATK's run_cpptraj does: at "\n", each line with trailing
// whitespace removed (str.rstrip), the rest kept for the next chunk; finish() gives the rest when
// it is not empty. Bytes that are not UTF-8 become U+FFFD (decode errors="replace"), also when a
// character is split between two chunks (MDATK decodes each 4 KiB chunk alone).
class LineSplitter {
public:
    std::vector<std::string> add(std::string_view chunk);
    std::optional<std::string> finish();

private:
    std::string m_pending;  // undecoded bytes of the current line
};

// " ".join(argv): how MDATK shows a command in its log and marker files.
std::string commandLineText(const std::vector<std::string>& argv);

// Python's subprocess.list2cmdline: the command line CreateProcessW gets on Windows.
std::string windowsCommandLine(const std::vector<std::string>& argv);

}  // namespace mdbridge::core::process
