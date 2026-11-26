// Running external programs: shutil.which, the child environment, a cancelable process runner
// with streamed output (see process.hpp).
#include "mdbridge/core/process.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <cwchar>
#include <system_error>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

#include "mdbridge/core/paths.hpp"

namespace mdbridge::core::process {

namespace fs = std::filesystem;

namespace {

#ifndef _WIN32
constexpr char kPathSeparator = ':';
#endif

std::vector<std::string> splitList(std::string_view text, char separator) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t end = text.find(separator, start);
        out.emplace_back(text.substr(start, end == std::string_view::npos ? std::string_view::npos : end - start));
        if (end == std::string_view::npos) {
            return out;
        }
        start = end + 1;
    }
}

#ifndef _WIN32
// path is `root` or inside it (both absolute, compared lexically)
bool inside(const std::string& path, const std::string& root) {
    if (path.empty() || root.empty()) {
        return false;
    }
    const fs::path p = pathFromUtf8(path).lexically_normal();
    const fs::path r = pathFromUtf8(root).lexically_normal();
    auto pi = p.begin();
    for (auto ri = r.begin(); ri != r.end(); ++ri) {
        if (ri->empty()) {
            continue;  // trailing separator
        }
        if (pi == p.end() || *pi != *ri) {
            return false;
        }
        ++pi;
    }
    return true;
}
#endif

// The bytes of `text` decoded as UTF-8 with errors="replace": each maximal invalid part (Unicode's
// "substitution of maximal subparts", which Python follows) becomes U+FFFD.
std::string replaceInvalidUtf8(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    const auto byte = [&text](std::size_t i) { return static_cast<unsigned char>(text[i]); };
    std::size_t i = 0;
    while (i < text.size()) {
        const unsigned char c = byte(i);
        if (c < 0x80) {
            out += static_cast<char>(c);
            ++i;
            continue;
        }
        std::size_t need = 0;
        unsigned char lo = 0x80, hi = 0xBF;  // allowed range of the second byte
        if (c >= 0xC2 && c <= 0xDF) {
            need = 1;
        } else if (c >= 0xE0 && c <= 0xEF) {
            need = 2;
            lo = c == 0xE0 ? 0xA0 : 0x80;
            hi = c == 0xED ? 0x9F : 0xBF;
        } else if (c >= 0xF0 && c <= 0xF4) {
            need = 3;
            lo = c == 0xF0 ? 0x90 : 0x80;
            hi = c == 0xF4 ? 0x8F : 0xBF;
        }
        if (need == 0) {
            out += "\xEF\xBF\xBD";
            ++i;
            continue;
        }
        std::size_t k = 1;
        while (k <= need && i + k < text.size()) {
            const unsigned char b = byte(i + k);
            const bool ok = k == 1 ? (b >= lo && b <= hi) : (b >= 0x80 && b <= 0xBF);
            if (!ok) {
                break;
            }
            ++k;
        }
        if (k == need + 1) {
            out.append(text.substr(i, k));
        } else {
            out += "\xEF\xBF\xBD";
        }
        i += k;
    }
    return out;
}

std::string rstrip(std::string text) {
    std::size_t end = text.size();
    while (end > 0) {
        const char c = text[end - 1];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f' || c == '\x1c' ||
            c == '\x1d' || c == '\x1e' || c == '\x1f') {
            --end;
        } else {
            break;
        }
    }
    text.resize(end);
    return text;
}

}  // namespace

// ---- environment --------------------------------------------------------------------------------

std::string environmentVariable(std::string_view name) {
#ifdef _WIN32
    const std::wstring wide = pathFromUtf8(name).wstring();
    const DWORD size = GetEnvironmentVariableW(wide.c_str(), nullptr, 0);
    if (size == 0) {
        return {};
    }
    std::wstring value(size, L'\0');
    const DWORD written = GetEnvironmentVariableW(wide.c_str(), value.data(), size);
    value.resize(written);
    return pathToUtf8(fs::path(value));
#else
    const char* value = std::getenv(std::string(name).c_str());
    return value ? value : "";
#endif
}

bool isExecutableFile(const std::string& path) {
    if (path.empty()) {
        return false;
    }
#ifdef _WIN32
    std::error_code ec;
    return fs::is_regular_file(pathFromUtf8(path), ec);
#else
    struct stat info {};
    return ::stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode) && ::access(path.c_str(), X_OK) == 0;
#endif
}

std::optional<std::string> which(std::string_view name) {
    const std::string cmd(name);
    if (cmd.empty()) {
        return std::nullopt;
    }
#ifdef _WIN32
    std::vector<std::string> extensions;
    {
        std::string pathext = environmentVariable("PATHEXT");
        if (pathext.empty()) {
            pathext = ".COM;.EXE;.BAT;.CMD;.VBS;.JS;.WS;.MSC";
        }
        for (const std::string& e : splitList(pathext, ';')) {
            if (!e.empty()) {
                extensions.push_back(e);
            }
        }
    }
    const auto lowerText = [](std::string t) {
        std::transform(t.begin(), t.end(), t.begin(),
                       [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; });
        return t;
    };
    std::vector<std::string> names;
    const std::string lowered = lowerText(cmd);
    const bool hasExtension = std::any_of(extensions.begin(), extensions.end(), [&](const std::string& e) {
        const std::string le = lowerText(e);
        return lowered.size() >= le.size() && lowered.compare(lowered.size() - le.size(), le.size(), le) == 0;
    });
    if (hasExtension) {
        names.push_back(cmd);
    } else {
        for (const std::string& e : extensions) {
            names.push_back(cmd + e);
        }
    }
    const auto check = [&names](const fs::path& folder) -> std::optional<std::string> {
        for (const std::string& n : names) {
            const fs::path candidate = folder.empty() ? pathFromUtf8(n) : folder / pathFromUtf8(n);
            if (isExecutableFile(pathToUtf8(candidate))) {
                return pathToUtf8(candidate);
            }
        }
        return std::nullopt;
    };
    if (cmd.find_first_of("/\\") != std::string::npos || (cmd.size() > 1 && cmd[1] == ':')) {
        return check({});
    }
    if (auto found = check(fs::current_path())) {
        return found;
    }
    for (const std::string& folder : splitList(environmentVariable("PATH"), ';')) {
        if (folder.empty()) {
            continue;
        }
        if (auto found = check(pathFromUtf8(folder))) {
            return found;
        }
    }
    return std::nullopt;
#else
    if (cmd.find('/') != std::string::npos) {
        return isExecutableFile(cmd) ? std::optional<std::string>(cmd) : std::nullopt;
    }
    const char* pathValue = std::getenv("PATH");
    const std::string path = pathValue ? pathValue : "/bin:/usr/bin";
    std::vector<std::string> seen;
    for (std::string folder : splitList(path, ':')) {
        if (folder.empty()) {
            folder = ".";  // os.curdir for an empty entry
        }
        if (std::find(seen.begin(), seen.end(), folder) != seen.end()) {
            continue;
        }
        seen.push_back(folder);
        const std::string candidate = folder.back() == '/' ? folder + cmd : folder + "/" + cmd;
        if (isExecutableFile(candidate)) {
            return candidate;
        }
    }
    return std::nullopt;
#endif
}

std::vector<std::pair<std::string, std::string>> childEnvironment() {
    std::vector<std::pair<std::string, std::string>> env;
#ifdef _WIN32
    LPWCH block = GetEnvironmentStringsW();
    if (block != nullptr) {
        for (const wchar_t* p = block; *p != L'\0'; p += std::wcslen(p) + 1) {
            const std::wstring entry(p);
            const std::size_t eq = entry.find(L'=', 1);  // "=C:=C:\\..." entries start with '='
            if (eq == std::wstring::npos) {
                continue;
            }
            env.emplace_back(pathToUtf8(fs::path(entry.substr(0, eq))), pathToUtf8(fs::path(entry.substr(eq + 1))));
        }
        FreeEnvironmentStringsW(block);
    }
#else
    for (char** p = environ; p != nullptr && *p != nullptr; ++p) {
        const std::string entry(*p);
        const std::size_t eq = entry.find('=');
        if (eq == std::string::npos || eq == 0) {
            continue;
        }
        env.emplace_back(entry.substr(0, eq), entry.substr(eq + 1));
    }
    std::string appdir;
    for (const auto& [name, value] : env) {
        if (name == "APPDIR") {
            appdir = value;
        }
    }
    if (!appdir.empty()) {  // started from the AppImage: drop what its start script added
        std::vector<std::pair<std::string, std::string>> kept;
        for (auto& [name, value] : env) {
            if (name == "FONTCONFIG_FILE" && inside(value, appdir)) {
                continue;
            }
            if (name == "LD_LIBRARY_PATH" || name == "PATH") {
                std::string joined;
                for (const std::string& part : splitList(value, kPathSeparator)) {
                    if (part.empty() || inside(part, appdir)) {
                        continue;
                    }
                    joined += (joined.empty() ? "" : std::string(1, kPathSeparator)) + part;
                }
                if (joined.empty()) {
                    continue;
                }
                value = joined;
            }
            kept.emplace_back(std::move(name), std::move(value));
        }
        env = std::move(kept);
    }
#endif
    return env;
}

// ---- texts ----------------------------------------------------------------------------------------

std::string commandLineText(const std::vector<std::string>& argv) {
    std::string out;
    for (std::size_t k = 0; k < argv.size(); ++k) {
        out += (k ? " " : "") + argv[k];
    }
    return out;
}

std::string windowsCommandLine(const std::vector<std::string>& argv) {
    std::string out;
    for (std::size_t k = 0; k < argv.size(); ++k) {
        const std::string& arg = argv[k];
        if (k) {
            out += ' ';
        }
        const bool quote = arg.empty() || arg.find_first_of(" \t") != std::string::npos;
        if (quote) {
            out += '"';
        }
        std::size_t backslashes = 0;
        for (char c : arg) {
            if (c == '\\') {
                ++backslashes;
            } else if (c == '"') {
                out.append(backslashes * 2 + 1, '\\');
                out += '"';
                backslashes = 0;
            } else {
                out.append(backslashes, '\\');
                backslashes = 0;
                out += c;
            }
        }
        out.append(backslashes * (quote ? 2 : 1), '\\');
        if (quote) {
            out += '"';
        }
    }
    return out;
}

std::vector<std::string> LineSplitter::add(std::string_view chunk) {
    m_pending.append(chunk);
    std::vector<std::string> lines;
    std::size_t start = 0;
    while (true) {
        const std::size_t nl = m_pending.find('\n', start);
        if (nl == std::string::npos) {
            break;
        }
        lines.push_back(rstrip(replaceInvalidUtf8(std::string_view(m_pending).substr(start, nl - start))));
        start = nl + 1;
    }
    m_pending.erase(0, start);
    return lines;
}

std::optional<std::string> LineSplitter::finish() {
    if (m_pending.empty()) {
        return std::nullopt;
    }
    std::string rest = rstrip(replaceInvalidUtf8(m_pending));
    m_pending.clear();
    return rest;
}

// ---- the runner ---------------------------------------------------------------------------------

#ifdef _WIN32

namespace {

std::string windowsErrorText(DWORD code, const std::string& what) {
    // Python's OSError texts for the common cases (FileNotFoundError: "[WinError 2] ...")
    const std::string message = std::system_category().message(static_cast<int>(code));
    return "[WinError " + std::to_string(code) + "] " + message + (what.empty() ? "" : ": '" + what + "'");
}

struct Handle {
    HANDLE h = nullptr;
    Handle() = default;
    explicit Handle(HANDLE value) : h(value) {}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    ~Handle() { reset(); }
    void reset(HANDLE value = nullptr) {
        if (h != nullptr && h != INVALID_HANDLE_VALUE) {
            CloseHandle(h);
        }
        h = value;
    }
};

}  // namespace

Result run(const std::vector<std::string>& argv, const RunOptions& options) {
    if (argv.empty()) {
        throw std::invalid_argument("process::run: empty command");
    }
    std::string program = argv[0];
    if (const std::optional<std::string> found = which(program)) {
        program = *found;
    }
    SECURITY_ATTRIBUTES inherit{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    HANDLE readEnd = nullptr, writeEnd = nullptr;
    if (!CreatePipe(&readEnd, &writeEnd, &inherit, 0)) {
        throw StartError(windowsErrorText(GetLastError(), ""));
    }
    Handle reader(readEnd), writer(writeEnd);
    SetHandleInformation(reader.h, HANDLE_FLAG_INHERIT, 0);
    Handle nul(CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &inherit, OPEN_EXISTING,
                           FILE_ATTRIBUTE_NORMAL, nullptr));

    STARTUPINFOW startup{};
    startup.cb = sizeof startup;
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = nul.h;
    startup.hStdOutput = writer.h;
    startup.hStdError = writer.h;

    // as Python's subprocess: no application name, the resolved program first on the command line
    // (CreateProcessW then also runs .bat / .cmd files, as conda's are on Windows)
    std::vector<std::string> args = argv;
    args[0] = program;
    std::wstring commandLine = pathFromUtf8(windowsCommandLine(args)).wstring();
    std::wstring block;
    const auto env = options.environment ? *options.environment : childEnvironment();
    for (const auto& [name, value] : env) {
        block += pathFromUtf8(name).wstring() + L"=" + pathFromUtf8(value).wstring();
        block += L'\0';
    }
    block += L'\0';
    const std::wstring cwd = options.workingDirectory.empty() ? std::wstring() : options.workingDirectory.wstring();

    Handle job(CreateJobObjectW(nullptr, nullptr));
    if (job.h != nullptr) {
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
        limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        SetInformationJobObject(job.h, JobObjectExtendedLimitInformation, &limits, sizeof limits);
    }
    PROCESS_INFORMATION info{};
    const BOOL started = CreateProcessW(nullptr, commandLine.data(), nullptr, nullptr, TRUE,
                                        CREATE_UNICODE_ENVIRONMENT | CREATE_NO_WINDOW | CREATE_SUSPENDED,
                                        block.data(), cwd.empty() ? nullptr : cwd.c_str(), &startup, &info);
    if (!started) {
        throw StartError(windowsErrorText(GetLastError(), program));
    }
    Handle process(info.hProcess), thread(info.hThread);
    if (job.h != nullptr) {
        AssignProcessToJobObject(job.h, process.h);
    }
    ResumeThread(thread.h);
    writer.reset();  // the child holds the write end now

    Result result;
    char buffer[4096];
    bool open = true;
    while (open) {
        if (!result.cancelled && options.cancelled && options.cancelled()) {
            result.cancelled = true;
            if (job.h != nullptr) {
                TerminateJobObject(job.h, 1);
            } else {
                TerminateProcess(process.h, 1);
            }
        }
        DWORD available = 0;
        if (!PeekNamedPipe(reader.h, nullptr, 0, nullptr, &available, nullptr)) {
            open = false;  // broken pipe: every writer has ended
            break;
        }
        if (available == 0) {
            if (WaitForSingleObject(process.h, 0) == WAIT_OBJECT_0) {
                // the program ended; grandchildren may still hold the pipe: keep reading until it breaks
                // unless nothing is left and the job is empty
                JOBOBJECT_BASIC_ACCOUNTING_INFORMATION accounting{};
                if (job.h == nullptr ||
                    !QueryInformationJobObject(job.h, JobObjectBasicAccountingInformation, &accounting,
                                               sizeof accounting, nullptr) ||
                    accounting.ActiveProcesses == 0) {
                    break;
                }
            }
            Sleep(static_cast<DWORD>(std::max<long long>(1, options.pollInterval.count())));
            continue;
        }
        DWORD got = 0;
        if (!ReadFile(reader.h, buffer, std::min<DWORD>(available, static_cast<DWORD>(sizeof buffer)), &got, nullptr) || got == 0) {
            break;
        }
        if (options.onOutput) {
            options.onOutput(std::string_view(buffer, got));
        }
    }
    WaitForSingleObject(process.h, INFINITE);
    DWORD code = 0;
    GetExitCodeProcess(process.h, &code);
    result.exitCode = static_cast<int>(code);
    return result;
}

#else

namespace {

std::string errnoText(int code, const std::string& what) {
    // Python's OSError str(): "[Errno 2] No such file or directory: 'cpptraj'"
    return "[Errno " + std::to_string(code) + "] " + std::strerror(code) + (what.empty() ? "" : ": '" + what + "'");
}

struct Fd {
    int fd = -1;
    Fd() = default;
    explicit Fd(int value) : fd(value) {}
    Fd(const Fd&) = delete;
    Fd& operator=(const Fd&) = delete;
    ~Fd() { reset(); }
    void reset(int value = -1) {
        if (fd >= 0) {
            ::close(fd);
        }
        fd = value;
    }
};

bool reap(pid_t pid, int& status, bool block) {
    while (true) {
        const pid_t r = ::waitpid(pid, &status, block ? 0 : WNOHANG);
        if (r == pid) {
            return true;
        }
        if (r == 0) {
            return false;
        }
        if (errno != EINTR) {
            status = 0;
            return true;
        }
    }
}

}  // namespace

Result run(const std::vector<std::string>& argv, const RunOptions& options) {
    if (argv.empty()) {
        throw std::invalid_argument("process::run: empty command");
    }
    std::string program = argv[0];
    if (program.find('/') == std::string::npos) {
        const std::optional<std::string> found = which(program);
        if (!found) {
            throw StartError(errnoText(ENOENT, program));
        }
        program = *found;
    }
    if (!options.workingDirectory.empty()) {
        std::error_code ec;
        if (!fs::is_directory(options.workingDirectory, ec)) {
            throw StartError(errnoText(ENOENT, pathToUtf8(options.workingDirectory)));
        }
    }
    // everything the child needs, prepared before fork (only async-signal-safe calls after it)
    std::vector<std::string> args = argv;
    std::vector<char*> argp;
    for (std::string& a : args) {
        argp.push_back(a.data());
    }
    argp.push_back(nullptr);
    std::vector<std::string> envText;
    for (const auto& [name, value] : options.environment ? *options.environment : childEnvironment()) {
        envText.push_back(name + "=" + value);
    }
    std::vector<char*> envp;
    for (std::string& e : envText) {
        envp.push_back(e.data());
    }
    envp.push_back(nullptr);
    const std::string cwd = options.workingDirectory.empty() ? std::string() : options.workingDirectory.string();

    int out[2] = {-1, -1};
    int err[2] = {-1, -1};
    if (::pipe2(out, O_CLOEXEC) != 0 || ::pipe2(err, O_CLOEXEC) != 0) {
        const int code = errno;
        for (int fd : {out[0], out[1], err[0], err[1]}) {
            if (fd >= 0) {
                ::close(fd);
            }
        }
        throw StartError(errnoText(code, ""));
    }
    Fd outRead(out[0]), outWrite(out[1]), errRead(err[0]), errWrite(err[1]);
    const pid_t pid = ::fork();
    if (pid < 0) {
        throw StartError(errnoText(errno, ""));
    }
    if (pid == 0) {
        ::setpgid(0, 0);  // its own process group: a cancel ends conda run's children too
        const int devnull = ::open("/dev/null", O_RDONLY);
        if (devnull >= 0) {
            ::dup2(devnull, 0);
        }
        ::dup2(out[1], 1);
        ::dup2(out[1], 2);
        int code = 0;
        if (!cwd.empty() && ::chdir(cwd.c_str()) != 0) {
            code = errno;
        } else {
            ::execve(program.c_str(), argp.data(), envp.data());
            code = errno;
        }
        const ssize_t ignored = ::write(err[1], &code, sizeof code);
        (void)ignored;
        ::_exit(127);
    }
    ::setpgid(pid, pid);  // also from the parent, whichever runs first
    outWrite.reset();
    errWrite.reset();
    int execError = 0;
    {
        ssize_t got = 0;
        do {
            got = ::read(errRead.fd, &execError, sizeof execError);
        } while (got < 0 && errno == EINTR);
        if (got == static_cast<ssize_t>(sizeof execError)) {
            int status = 0;
            reap(pid, status, true);
            throw StartError(errnoText(execError, program));
        }
    }

    Result result;
    int status = 0;
    bool exited = false;
    std::chrono::steady_clock::time_point killAt{};
    bool killed = false;
    char buffer[4096];
    const int timeout = static_cast<int>(std::max<long long>(1, options.pollInterval.count()));
    while (true) {
        if (!result.cancelled && options.cancelled && options.cancelled()) {
            result.cancelled = true;
            ::kill(-pid, SIGTERM);
            killAt = std::chrono::steady_clock::now() + options.terminateGrace;
        }
        if (result.cancelled && !killed && std::chrono::steady_clock::now() >= killAt) {
            ::kill(-pid, SIGKILL);
            killed = true;
        }
        pollfd p{outRead.fd, POLLIN, 0};
        const int ready = ::poll(&p, 1, timeout);
        if (ready < 0) {
            if (errno == EINTR) {
                continue;
            }
            break;
        }
        if (ready == 0) {
            if (!exited) {
                exited = reap(pid, status, false);
            }
            if (exited && result.cancelled) {
                break;  // a cancelled program's leftover children may keep the pipe open
            }
            continue;
        }
        const ssize_t got = ::read(outRead.fd, buffer, sizeof buffer);
        if (got < 0) {
            if (errno == EINTR || errno == EAGAIN) {
                continue;
            }
            break;
        }
        if (got == 0) {
            break;  // every writer has closed the pipe
        }
        if (options.onOutput) {
            options.onOutput(std::string_view(buffer, static_cast<std::size_t>(got)));
        }
    }
    if (!exited) {
        if (result.cancelled && !killed) {
            // wait for the grace period, then kill
            while (!(exited = reap(pid, status, false)) && std::chrono::steady_clock::now() < killAt) {
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            }
            if (!exited) {
                ::kill(-pid, SIGKILL);
            }
        }
        if (!exited) {
            reap(pid, status, true);
        }
    }
    if (result.cancelled) {
        ::kill(-pid, SIGKILL);  // whatever of its group is left
    }
    result.exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : WIFSIGNALED(status) ? -WTERMSIG(status) : 0;
    return result;
}

#endif

}  // namespace mdbridge::core::process
