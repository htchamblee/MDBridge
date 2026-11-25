// Conda environments that provide cpptraj (mdatk/cleaning/conda.py candidate_roots,
// _registered_prefixes, find_envs, find_cpptraj_envs).
#include "mdbridge/core/conda.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <system_error>
#include <tuple>

#ifndef _WIN32
#include <unistd.h>
#endif

#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#ifdef _WIN32
#include "mdbridge/core/ntpath.hpp"
#endif

namespace mdbridge::core::conda {

namespace fs = std::filesystem;

namespace {

// os.environ.get(name) ("" when unset)
std::string environment(const char* name) {
#ifdef _WIN32
    std::size_t size = 0;
    if (_wgetenv_s(&size, nullptr, 0, fs::path(name).c_str()) != 0 || size == 0) {
        return {};
    }
    std::wstring value(size, L'\0');
    if (_wgetenv_s(&size, value.data(), size, fs::path(name).c_str()) != 0) {
        return {};
    }
    value.resize(size ? size - 1 : 0);
    return pathToUtf8(fs::path(value));
#else
    const char* value = std::getenv(name);
    return value ? value : "";
#endif
}

// os.path of the platform
std::string expanduser(std::string_view path) {
#ifdef _WIN32
    return ntpath::expanduser(path);
#else
    return posixpath::expanduser(path);
#endif
}

std::string dirname(std::string_view path) {
#ifdef _WIN32
    return ntpath::dirname(path);
#else
    return posixpath::dirname(path);
#endif
}

std::string basename(std::string_view path) {
#ifdef _WIN32
    return pathBasename(path);
#else
    return posixpath::basename(path);
#endif
}

std::string join(std::string_view base, std::string_view name) {
#ifdef _WIN32
    return pathToUtf8(pathFromUtf8(base) / pathFromUtf8(name));
#else
    return posixpath::join(base, name);
#endif
}

// os.path.realpath (non-strict): absolute, normalised, symbolic links of existing parts resolved
std::string realpath(std::string_view path) {
#ifdef _WIN32
    std::error_code ec;
    fs::path absolute = fs::absolute(pathFromUtf8(path), ec);
    if (ec) {
        absolute = pathFromUtf8(ntpath::abspath(path, currentDirectory()));
    }
    fs::path real = fs::weakly_canonical(absolute, ec);
    if (ec) {
        real = absolute.lexically_normal();
    }
    return pathToUtf8(real.make_preferred());
#else
    const std::string absolute = posixpath::abspath(path, currentDirectory());
    std::error_code ec;
    const fs::path real = fs::weakly_canonical(pathFromUtf8(absolute), ec);
    if (ec) {
        return absolute;
    }
    std::string out = pathToUtf8(real);
    while (out.size() > 1 && out.back() == '/') {  // weakly_canonical keeps a trailing separator
        out.pop_back();
    }
    return out;
#endif
}

bool isDirectory(const std::string& path) {
    std::error_code ec;
    return fs::is_directory(pathFromUtf8(path), ec);
}

// The root of an environment prefix: <root>/envs/<name> -> <root>, else the prefix itself.
std::string rootOf(const std::string& prefix) {
    const std::string parent = dirname(prefix);
    return basename(parent) == "envs" ? dirname(parent) : prefix;
}

// _registered_prefixes: the non-blank lines of ~/.conda/environments.txt, stripped
std::vector<std::string> registeredPrefixes() {
    std::ifstream in(pathFromUtf8(expanduser("~/.conda/environments.txt")), std::ios::binary);
    if (!in) {
        return {};
    }
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {  // universal newlines: \n, \r\n and \r end a line
        std::size_t end = text.find_first_of("\r\n", start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string line = pythonStrip(std::string_view(text).substr(start, end - start));
        if (!line.empty()) {
            out.push_back(std::move(line));
        }
        start = end + 1;
        if (end + 1 < text.size() && text[end] == '\r' && text[end + 1] == '\n') {
            ++start;
        }
    }
    return out;
}

bool containsAmber(const std::string& name) {  // re.search("amber", name, re.IGNORECASE)
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return lower.find("amber") != std::string::npos;
}

}  // namespace

std::string CondaEnv::cpptraj() const { return join(join(prefix, "bin"), "cpptraj"); }

bool CondaEnv::hasCpptraj() const {
    const std::string path = cpptraj();
    std::error_code ec;
    if (!fs::is_regular_file(pathFromUtf8(path), ec)) {
        return false;
    }
#ifdef _WIN32
    return true;  // os.access(path, X_OK) on Windows: the file exists
#else
    return ::access(path.c_str(), X_OK) == 0;
#endif
}

std::vector<std::string> candidateRoots(const std::vector<std::string>& extra) {
    std::vector<std::string> roots(extra);
    if (const std::string condaExe = environment("CONDA_EXE"); !condaExe.empty()) {
        roots.push_back(dirname(dirname(condaExe)));
    }
    if (const std::string prefix = environment("CONDA_PREFIX"); !prefix.empty()) {
        roots.push_back(rootOf(prefix));
    }
    for (const char* variable : {"MAMBA_ROOT_PREFIX", "CONDA_ROOT"}) {
        if (std::string value = environment(variable); !value.empty()) {
            roots.push_back(std::move(value));
        }
    }
    for (std::string_view root : kCommonRoots) {
        roots.emplace_back(root);
    }
    std::vector<std::string> out;
    std::set<std::string> seen;
    for (const std::string& r : roots) {
        std::string path = realpath(expanduser(r));
        if (seen.count(path) == 0 && isDirectory(path)) {
            seen.insert(path);
            out.push_back(std::move(path));
        }
    }
    return out;
}

std::vector<CondaEnv> findEnvs(const std::optional<std::vector<std::string>>& givenRoots, bool includeRegistered) {
    std::vector<std::string> roots;
    if (givenRoots) {
        for (const std::string& r : *givenRoots) {
            roots.push_back(realpath(expanduser(r)));
        }
    } else {
        roots = candidateRoots();
    }
    std::vector<CondaEnv> envs;
    std::set<std::string> seen;
    const auto add = [&](const std::string& name, const std::string& prefix, const std::string& root) {
        std::string real = realpath(prefix);
        if (seen.count(real) == 0 && isDirectory(real)) {
            seen.insert(real);
            envs.push_back(CondaEnv{name, std::move(real), root});
        }
    };
    for (const std::string& root : roots) {
        add("base", root, root);
        const std::string envDir = join(root, "envs");
        if (isDirectory(envDir)) {
            std::vector<std::string> names;  // os.listdir, sorted
            std::error_code ec;
            for (fs::directory_iterator it(pathFromUtf8(envDir), ec), end; !ec && it != end; it.increment(ec)) {
                names.push_back(pathToUtf8(it->path().filename()));
            }
            std::sort(names.begin(), names.end());
            for (const std::string& name : names) {
                add(name, join(envDir, name), root);
            }
        }
    }
    if (includeRegistered) {
        for (const std::string& prefix : registeredPrefixes()) {
            const std::string root = rootOf(prefix);
            add(prefix != root ? basename(prefix) : "base", prefix, root);
        }
    }
    return envs;
}

std::vector<CondaEnv> findCpptrajEnvs(const std::optional<std::vector<std::string>>& roots) {
    std::vector<CondaEnv> envs;
    for (CondaEnv& env : findEnvs(roots, !roots.has_value())) {
        if (env.hasCpptraj()) {
            envs.push_back(std::move(env));
        }
    }
    std::stable_sort(envs.begin(), envs.end(), [](const CondaEnv& a, const CondaEnv& b) {
        return std::make_tuple(!containsAmber(a.name), a.name == "base", std::string_view(a.name)) <
               std::make_tuple(!containsAmber(b.name), b.name == "base", std::string_view(b.name));
    });
    return envs;
}

}  // namespace mdbridge::core::conda
