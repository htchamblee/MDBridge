// Conda environments with cpptraj (mdatk/cleaning/conda.py find_cpptraj_envs) on fake
// installations: environments recognised by an executable bin/cpptraj, names with "amber" first,
// "base" last, real paths listed once, roots from the environment variables MDATK reads and from
// ~/.conda/environments.txt. tests/parity/cli_parity.py (`list_conda`) compares the same discovery
// with MDATK's on Linux.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/conda.hpp"
#include "mdbridge/core/paths.hpp"
#include "test_support.hpp"

namespace mdbridge {
namespace {

namespace fs = std::filesystem;
using core::conda::CondaEnv;

// Environment variables set for a while.
class ScopedEnvironment {
public:
    ScopedEnvironment() = default;
    ScopedEnvironment(const ScopedEnvironment&) = delete;
    ScopedEnvironment& operator=(const ScopedEnvironment&) = delete;
    ~ScopedEnvironment() {
        for (auto it = m_saved.rbegin(); it != m_saved.rend(); ++it) {
            put(it->first, it->second);
        }
    }
    void set(const std::string& name, const std::string& value) {
        m_saved.emplace_back(name, get(name));
        put(name, value);
    }

private:
    static std::optional<std::string> get(const std::string& name) {
#ifdef _WIN32
        char* value = nullptr;
        std::size_t size = 0;
        if (_dupenv_s(&value, &size, name.c_str()) != 0 || value == nullptr) {
            return std::nullopt;
        }
        std::string out(value);
        std::free(value);
        return out;
#else
        const char* value = std::getenv(name.c_str());
        return value ? std::optional<std::string>(value) : std::nullopt;
#endif
    }
    static void put(const std::string& name, const std::optional<std::string>& value) {
#ifdef _WIN32
        _putenv_s(name.c_str(), value ? value->c_str() : "");
#else
        if (value) {
            setenv(name.c_str(), value->c_str(), 1);
        } else {
            unsetenv(name.c_str());
        }
#endif
    }
    std::vector<std::pair<std::string, std::optional<std::string>>> m_saved;
};

// os.path.realpath as the platform returns it ('\' separators on Windows)
std::string real(const fs::path& path) {
    fs::path p = fs::canonical(path);
    return core::pathToUtf8(p.make_preferred());
}

void makeExecutable(const fs::path& path) {
    fs::create_directories(path.parent_path());
    test::writeFile(path, "#!/bin/sh\n");
    fs::permissions(path, fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec,
                    fs::perm_options::add);
}

std::vector<std::pair<std::string, std::string>> namesAndPrefixes(const std::vector<CondaEnv>& envs,
                                                                   const fs::path& under) {
    std::vector<std::pair<std::string, std::string>> out;
    const std::string base = real(under);
    for (const CondaEnv& e : envs) {
        if (e.prefix.starts_with(base)) {  // installations of the machine itself are left out
            out.emplace_back(e.name, e.prefix);
        }
    }
    return out;
}

TEST(Conda, EnvironmentsWithCpptrajUnderGivenRoots) {
    test::TempDir work;
    const fs::path root = work.path() / "miniforge";
    makeExecutable(root / "bin" / "cpptraj");
    makeExecutable(root / "envs" / "zeta" / "bin" / "cpptraj");
    makeExecutable(root / "envs" / "AmberTools26" / "bin" / "cpptraj");
    makeExecutable(root / "envs" / "my-amber" / "bin" / "cpptraj");
    fs::create_directories(root / "envs" / "py312" / "bin");           // no cpptraj
    fs::create_directories(root / "envs" / "folder" / "bin" / "cpptraj");  // not a file
    test::writeFile(root / "envs" / "stray-file", "");                    // not an environment
#ifndef _WIN32
    fs::create_directories(root / "envs" / "noexec" / "bin");
    test::writeFile(root / "envs" / "noexec" / "bin" / "cpptraj", "");  // not executable
    fs::create_directory_symlink(root / "envs" / "zeta", root / "envs" / "zeta-link");  // listed once
#endif
    const std::vector<CondaEnv> envs = core::conda::findCpptrajEnvs(std::vector<std::string>{
        core::pathToUtf8(root), core::pathToUtf8(root), core::pathToUtf8(work.path() / "none")});
    const std::vector<std::pair<std::string, std::string>> expected{
        {"AmberTools26", real(root / "envs" / "AmberTools26")},
        {"my-amber", real(root / "envs" / "my-amber")},
        {"zeta", real(root / "envs" / "zeta")},
        {"base", real(root)},
    };
    EXPECT_EQ(namesAndPrefixes(envs, work.path()), expected);
    for (const CondaEnv& e : envs) {
        EXPECT_EQ(e.root, real(root));
        EXPECT_TRUE(e.hasCpptraj());
        EXPECT_EQ(core::pathFromUtf8(e.cpptraj()), core::pathFromUtf8(e.prefix) / "bin" / "cpptraj");
    }

    // every environment (find_envs), without the registered ones
    std::vector<std::string> names;
    for (const CondaEnv& e : core::conda::findEnvs(std::vector<std::string>{core::pathToUtf8(root)}, false)) {
        names.push_back(e.name);
    }
#ifdef _WIN32
    EXPECT_EQ(names, (std::vector<std::string>{"base", "AmberTools26", "folder", "my-amber", "py312", "zeta"}));
#else
    EXPECT_EQ(names,
              (std::vector<std::string>{"base", "AmberTools26", "folder", "my-amber", "noexec", "py312", "zeta"}));
#endif
}

TEST(Conda, RootsFromTheEnvironmentAndRegisteredEnvironments) {
    test::TempDir work;
    const fs::path home = work.path() / "home";
    const fs::path forge = work.path() / "miniforge";
    const fs::path mamba = work.path() / "mamba";
    const fs::path other = work.path() / "other";
    makeExecutable(forge / "envs" / "amber24" / "bin" / "cpptraj");
    makeExecutable(mamba / "envs" / "tools" / "bin" / "cpptraj");
    makeExecutable(work.path() / "croot" / "bin" / "cpptraj");
    makeExecutable(other / "envs" / "AmberTools25" / "bin" / "cpptraj");
    makeExecutable(work.path() / "standalone" / "bin" / "cpptraj");
    fs::create_directories(forge / "condabin");
    fs::create_directories(home / ".conda");
    test::writeFile(home / ".conda" / "environments.txt",
                    "  " + core::pathToUtf8(other / "envs" / "AmberTools25") + "\r\n\n" +
                        core::pathToUtf8(work.path() / "standalone") + "\r" +
                        core::pathToUtf8(work.path() / "missing") + "\n" +
                        core::pathToUtf8(forge / "envs" / "amber24") + "\n");

    ScopedEnvironment env;
#ifdef _WIN32
    env.set("USERPROFILE", core::pathToUtf8(home));
#else
    env.set("HOME", core::pathToUtf8(home));
#endif
    env.set("CONDA_EXE", core::pathToUtf8(forge / "condabin" / "conda"));
    env.set("CONDA_PREFIX", core::pathToUtf8(mamba / "envs" / "tools"));  // <root>/envs/<name>: the root
    env.set("MAMBA_ROOT_PREFIX", "");                                     // empty: unset
    env.set("CONDA_ROOT", core::pathToUtf8(work.path() / "croot"));

    const std::vector<std::string> roots = core::conda::candidateRoots({core::pathToUtf8(forge)});
    std::vector<std::string> mine;
    for (const std::string& r : roots) {
        if (r.starts_with(real(work.path()))) {
            mine.push_back(r);
        }
    }
    // extra first, then CONDA_EXE's root (the same: once), CONDA_PREFIX's, CONDA_ROOT, ~/.conda
    EXPECT_EQ(mine, (std::vector<std::string>{real(forge), real(mamba), real(work.path() / "croot"),
                                              real(home / ".conda")}));

    const std::vector<std::pair<std::string, std::string>> expected{
        {"AmberTools25", real(other / "envs" / "AmberTools25")},
        {"amber24", real(forge / "envs" / "amber24")},
        {"tools", real(mamba / "envs" / "tools")},
        {"base", real(work.path() / "croot")},
        {"base", real(work.path() / "standalone")},
    };
    EXPECT_EQ(namesAndPrefixes(core::conda::findCpptrajEnvs(), work.path()), expected);
}

}  // namespace
}  // namespace mdbridge
