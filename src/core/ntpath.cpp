#include "mdbridge/core/ntpath.hpp"

#include <cstdlib>
#include <vector>

#ifdef _WIN32
#include <filesystem>

#include "mdbridge/core/paths.hpp"
#endif

namespace mdbridge::core::ntpath {

namespace {

constexpr bool isAsciiLetter(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }

constexpr char asciiLower(char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; }

// "C:" (an ASCII letter and a colon)
bool isLetterDrive(std::string_view drive) {
    return drive.size() == 2 && isAsciiLetter(drive[0]) && drive[1] == ':';
}

// The same drive, letters compared without case (UNC drives as given)
bool sameDrive(std::string_view a, std::string_view b) {
    if (isLetterDrive(a) && isLetterDrive(b)) {
        return asciiLower(a[0]) == asciiLower(b[0]);
    }
    return a == b;
}

std::size_t findSep(std::string_view path, std::size_t from) {
    for (std::size_t i = from; i < path.size(); ++i) {
        if (isSep(path[i])) {
            return i;
        }
    }
    return std::string_view::npos;
}

#ifdef _WIN32
std::optional<std::string> environment(const wchar_t* name) {
    std::size_t size = 0;
    if (_wgetenv_s(&size, nullptr, 0, name) != 0 || size == 0) {
        return std::nullopt;
    }
    std::wstring value(size, L'\0');
    if (_wgetenv_s(&size, value.data(), value.size(), name) != 0) {
        return std::nullopt;
    }
    value.resize(size ? size - 1 : 0);
    return pathToUtf8(std::filesystem::path(value));
}
#else
std::optional<std::string> environment(const char* name) {
    const char* value = std::getenv(name);
    return value ? std::optional<std::string>(value) : std::nullopt;
}
#endif

}  // namespace

Parts splitroot(std::string_view path) {
    if (!path.empty() && isSep(path[0])) {
        if (path.size() >= 2 && isSep(path[1])) {
            // UNC: "\\server\share" (or "\\?\UNC\server\share"), then the root
            std::size_t start = 2;
            if (path.size() >= 8) {
                const std::string_view head = path.substr(0, 8);
                const bool unc = isSep(head[0]) && isSep(head[1]) && head[2] == '?' && isSep(head[3]) &&
                                 (head[4] == 'U' || head[4] == 'u') && (head[5] == 'N' || head[5] == 'n') &&
                                 (head[6] == 'C' || head[6] == 'c') && isSep(head[7]);
                if (unc) {
                    start = 8;
                }
            }
            const std::size_t index = findSep(path, start);
            if (index == std::string_view::npos) {
                return {path, {}, {}};
            }
            const std::size_t index2 = findSep(path, index + 1);
            if (index2 == std::string_view::npos) {
                return {path, {}, {}};
            }
            return {path.substr(0, index2), path.substr(index2, 1), path.substr(index2 + 1)};
        }
        return {{}, path.substr(0, 1), path.substr(1)};
    }
    if (path.size() >= 2 && isAsciiLetter(path[0]) && path[1] == ':') {
        if (path.size() >= 3 && isSep(path[2])) {
            return {path.substr(0, 2), path.substr(2, 1), path.substr(3)};
        }
        return {path.substr(0, 2), {}, path.substr(2)};
    }
    return {{}, {}, path};
}

bool isDevicePath(std::string_view path) {
    return path.size() >= 4 && isSep(path[0]) && isSep(path[1]) && (path[2] == '?' || path[2] == '.') &&
           isSep(path[3]);
}

bool isabs(std::string_view path) {
    const Parts parts = splitroot(path);
    return !parts.root.empty() || (!parts.drive.empty() && !isLetterDrive(parts.drive));
}

std::string generic(std::string_view path) {
    std::string out(path);
    for (char& c : out) {
        if (c == '\\') {
            c = '/';
        }
    }
    return out;
}

std::string join(std::string_view base, std::string_view path) {
    const Parts b = splitroot(base);
    const Parts p = splitroot(path);
    std::string drive(b.drive);
    std::string root(b.root);
    std::string tail(b.tail);
    if (!p.root.empty()) {
        if (!p.drive.empty() || drive.empty()) {
            drive = p.drive;
        }
        root = p.root;
        tail = p.tail;
    } else if (!p.drive.empty() && !sameDrive(p.drive, drive)) {
        return std::string(path);
    } else {
        if (!p.drive.empty()) {
            drive = p.drive;  // the same drive, as the path writes it
        }
        if (!tail.empty() && !isSep(tail.back())) {
            tail += '\\';
        }
        tail += p.tail;
    }
    // a UNC drive followed by a tail needs a separator between them
    if (!tail.empty() && root.empty() && !drive.empty() && drive.back() != ':' && !isSep(drive.back())) {
        return drive + '\\' + tail;
    }
    return drive + root + tail;
}

std::string normpath(std::string_view path) {
    if (isDevicePath(path)) {
        return std::string(path);
    }
    const Parts parts = splitroot(path);
    std::string out = generic(parts.drive);
    if (!parts.root.empty()) {
        out += '/';
    }
    const bool rooted = !parts.root.empty();
    std::vector<std::string_view> comps;
    std::size_t start = 0;
    const std::string_view tail = parts.tail;
    while (start <= tail.size()) {
        std::size_t end = findSep(tail, start);
        if (end == std::string_view::npos) {
            end = tail.size();
        }
        const std::string_view comp = tail.substr(start, end - start);
        start = end + 1;
        if (comp.empty() || comp == ".") {
            continue;
        }
        if (comp == "..") {
            if (!comps.empty() && comps.back() != "..") {
                comps.pop_back();
                continue;
            }
            if (comps.empty() && rooted) {
                continue;  // never above the root
            }
        }
        comps.push_back(comp);
    }
    if (out.empty() && comps.empty()) {
        return ".";
    }
    for (std::size_t k = 0; k < comps.size(); ++k) {
        if (k) {
            out += '/';
        }
        out += comps[k];
    }
    return out;
}

std::string abspath(std::string_view path, std::string_view base) {
    if (isDevicePath(path)) {
        return std::string(path);
    }
    std::string joined = join(base, path);
    const Parts parts = splitroot(joined);
    if (isLetterDrive(parts.drive) && parts.root.empty()) {  // "C:x" not on the base's drive
        joined.insert(2, 1, '/');
    }
    return normpath(joined);
}

std::string dirname(std::string_view path) {
    const Parts parts = splitroot(path);
    std::size_t i = parts.tail.size();
    while (i > 0 && !isSep(parts.tail[i - 1])) {
        --i;
    }
    while (i > 0 && isSep(parts.tail[i - 1])) {
        --i;
    }
    std::string out(parts.drive);
    out += parts.root;
    out += parts.tail.substr(0, i);
    return out;
}

std::optional<std::string> profileHome(const std::optional<std::string>& userProfile,
                                       const std::optional<std::string>& homeDrive,
                                       const std::optional<std::string>& homePath) {
    if (userProfile) {
        return userProfile;
    }
    if (!homePath) {
        return std::nullopt;
    }
    return join(homeDrive.value_or(std::string()), *homePath);
}

std::optional<std::string> environmentHome() {
#ifdef _WIN32
    return profileHome(environment(L"USERPROFILE"), environment(L"HOMEDRIVE"), environment(L"HOMEPATH"));
#else
    return profileHome(environment("USERPROFILE"), environment("HOMEDRIVE"), environment("HOMEPATH"));
#endif
}

std::string expanduser(std::string_view path, const std::optional<std::string>& home) {
    if (path.empty() || path.front() != '~' || (path.size() > 1 && !isSep(path[1])) || !home) {
        return std::string(path);
    }
    std::string out = *home;
    out += path.substr(1);
    return out;
}

std::string expanduser(std::string_view path) {
    if (path.empty() || path.front() != '~') {
        return std::string(path);
    }
    return expanduser(path, environmentHome());
}

std::string resolve(std::string_view path, std::string_view baseDir, const std::optional<std::string>& home) {
    if (path.empty()) {
        return {};
    }
    return abspath(expanduser(path, home), baseDir);
}

std::string resolve(std::string_view path, std::string_view baseDir) {
    if (path.empty()) {
        return {};
    }
    return abspath(expanduser(path), baseDir);
}

}  // namespace mdbridge::core::ntpath
