#pragma once
// Conda environments that provide cpptraj (mdatk/cleaning/conda.py: candidate_roots, find_envs,
// find_cpptraj_envs). People install AmberTools under different environment names
// (AmberTools26, ambertools, amber24 ...), so environments are recognised by what they contain - an
// executable <prefix>/bin/cpptraj - rather than by their name.
//
// Paths are UTF-8 and follow the platform's os.path: on Linux and macOS posixpath ("~" from $HOME,
// else the passwd entry), on Windows ntpath ("~" from %USERPROFILE%, else %HOMEDRIVE%%HOMEPATH%;
// results with '\' separators, as os.path.realpath returns them there; every existing file counts
// as executable, as os.access(path, X_OK) does on Windows). Real paths resolve symbolic links of
// the parts that exist (os.path.realpath, non-strict).

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::core::conda {

// COMMON_ROOTS: installation folders searched besides those the environment variables name.
inline constexpr std::string_view kCommonRoots[] = {
    "~/miniforge3",       "~/mambaforge",          "~/miniconda3",          "~/anaconda3",
    "~/micromamba",       "~/.local/share/mamba",  "~/.conda",              "/opt/conda",
    "/opt/miniforge3",    "/opt/mambaforge",       "/opt/miniconda3",       "/opt/anaconda3",
    "/usr/local/miniconda3", "/usr/local/anaconda3",
};

struct CondaEnv {
    std::string name;    // environment name ("base" for a root installation)
    std::string prefix;  // absolute environment directory (real path)
    std::string root;    // installation the environment belongs to

    // <prefix>/bin/cpptraj
    std::string cpptraj() const;
    // cpptraj() is a file and executable
    bool hasCpptraj() const;
};

// candidate_roots(extra): `extra`, then the root of $CONDA_EXE (<root>/bin/conda or
// <root>/condabin/conda), of $CONDA_PREFIX (<root>/envs/<name>, or a root itself),
// $MAMBA_ROOT_PREFIX, $CONDA_ROOT and kCommonRoots; "~" expanded, real paths, existing folders
// only, each once. Empty variables count as unset.
std::vector<std::string> candidateRoots(const std::vector<std::string>& extra = {});

// find_envs(roots, include_registered): every root ("base") and the folders of its envs/ (sorted
// by name), then, with includeRegistered, the environments listed in ~/.conda/environments.txt;
// existing folders only, each real path once. `roots`: candidateRoots() when not given.
std::vector<CondaEnv> findEnvs(const std::optional<std::vector<std::string>>& roots = std::nullopt,
                               bool includeRegistered = true);

// find_cpptraj_envs(roots): the environments with an executable cpptraj (registered ones too when
// `roots` is not given), names containing "amber" (any case) first, then the others, "base" last
// within each part, by name within that (stable).
std::vector<CondaEnv> findCpptrajEnvs(const std::optional<std::vector<std::string>>& roots = std::nullopt);

// detect_ambertools_env(roots): the first of findCpptrajEnvs(roots).
std::optional<CondaEnv> detectAmbertoolsEnv(const std::optional<std::vector<std::string>>& roots = std::nullopt);

// resolve_env(spec, roots): an environment given by a path (spec contains the platform's
// separator: its real path when it is a folder; the root is two levels up for <root>/envs/<name>,
// else the prefix itself) or by name (the first of findEnvs(roots) with that name). Blank: none.
std::optional<CondaEnv> resolveEnv(std::string_view spec,
                                   const std::optional<std::vector<std::string>>& roots = std::nullopt);

// find_conda_exe(explicit, env): `explicit` when given (a path: when executable; a name: found on
// PATH); otherwise <root>/condabin/conda, bin/conda, bin/mamba, bin/micromamba of env's root and
// $CONDA_EXE, the first executable file, then conda, mamba, micromamba on PATH. "~" expanded.
std::optional<std::string> findCondaExe(std::string_view explicitExe = {}, const std::optional<CondaEnv>& env = std::nullopt);

// run_command(conda_exe, env_spec, argv): [exe, "run", ("--no-capture-output",) "-p"|"-n", env, argv...]
// ("-p" when env_spec contains the platform's separator; micromamba gets no --no-capture-output).
std::vector<std::string> runCommand(const std::string& condaExe, const std::string& envSpec,
                                    const std::vector<std::string>& argv);

}  // namespace mdbridge::core::conda
