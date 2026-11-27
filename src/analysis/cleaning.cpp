// Trajectory cleaning with cpptraj: a port of MDATK's mdatk/cleaning/cpptraj.py (MDATK's own code).
#include "mdbridge/analysis/cleaning.hpp"

#include <algorithm>
#include <cstdio>
#include <optional>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <system_error>

#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/conda.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/process.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/version.hpp"
#ifdef _WIN32
#include "mdbridge/core/ntpath.hpp"
#else
#include <unistd.h>
#endif

namespace mdbridge::analysis {

namespace fs = std::filesystem;
using core::Table;
using core::Value;

namespace {

std::string abspath(const std::string& path) {
#ifdef _WIN32
    return core::ntpath::abspath(path, core::currentDirectory());
#else
    return core::posixpath::abspath(path, core::currentDirectory());
#endif
}

std::string joinPath(const std::string& dir, std::string_view name) {
#ifdef _WIN32
    return core::pathToUtf8(core::pathFromUtf8(dir) / core::pathFromUtf8(name));
#else
    return core::posixpath::join(dir, name);
#endif
}

std::string expanduser(const std::string& path) {
#ifdef _WIN32
    return core::ntpath::expanduser(path);
#else
    return core::posixpath::expanduser(path);
#endif
}

// os.sep in path; on Windows '/' too (DIFFERENCES "Windows paths" 4: cpptraj_path with '/' is a path)
bool hasSeparator(std::string_view text) {
#ifdef _WIN32
    return text.find_first_of("/\\") != std::string_view::npos;
#else
    return text.find('/') != std::string_view::npos;
#endif
}

bool exists(const std::string& path) {
    std::error_code ec;
    return !path.empty() && fs::exists(core::pathFromUtf8(path), ec);
}

bool isParm7(const std::string& topology) {
    std::string ext = core::posixpath::splitext(topology).second;
    for (char& c : ext) {
        c = (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
    }
    return ext == ".parm7" || ext == ".prmtop" || ext == ".top";
}

// datetime.now().isoformat(timespec="seconds")
std::string isoNow() {
    const core::DateTime now = core::clockNow();
    char text[64];
    std::snprintf(text, sizeof text, "%04d-%02d-%02dT%02d:%02d:%02d", now.year, now.month, now.day, now.hour,
                  now.minute, now.second);
    return text;
}

std::string quote(const std::string& path) {
    if (path.find('"') != std::string::npos) {
        throw CpptrajError("path contains a double quote: " + path);
    }
    return "\"" + path + "\"";
}

// _compact_mask: cpptraj splits command arguments on whitespace; masks never need it
std::string compactMask(const std::string& mask) {
    std::u32string out;
    for (const char32_t c : core::decodeUtf8Surrogateescape(mask)) {
        if (!core::isPythonSpace(c)) {
            out += c;
        }
    }
    return core::encodeUtf8(out);
}

std::string stripMask(const core::CleaningParams& params) {
    std::vector<std::string> masks;
    for (const std::string* m : {&params.solvent_mask, &params.ion_mask, &params.extra_strip_mask}) {
        std::string c = compactMask(*m);
        if (!c.empty()) {
            masks.push_back(std::move(c));
        }
    }
    if (masks.size() > 1) {
        std::string out;
        for (std::size_t k = 0; k < masks.size(); ++k) {
            out += (k ? "|(" : "(") + masks[k] + ")";
        }
        return out;
    }
    return masks.empty() ? std::string() : masks.front();
}

// _PERCENT.findall(text)[-1]: the last (?<![\d.])(\d{1,3})% of the text
std::optional<int> lastPercent(std::string_view text) {
    std::optional<int> last;
    const auto digit = [](char c) { return c >= '0' && c <= '9'; };
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] != '%') {
            continue;
        }
        std::size_t j = i;
        while (j > 0 && digit(text[j - 1])) {
            --j;
        }
        const std::size_t length = i - j;
        if (length >= 1 && length <= 3 && (j == 0 || text[j - 1] != '.')) {
            last = std::stoi(std::string(text.substr(j, length)));
        }
    }
    return last;
}

// The cleaning settings as vars(params) gives them (every field, in the model's order).
Value paramsTable(const core::CleaningParams& params) {
    core::Project holder;
    holder.cleaning = params;
    const Value plain = core::toPlain(holder);
    const Value* section = plain.find("cleaning");
    return section ? *section : Value(Table{});
}

void writeTextFile(const std::string& path, const std::string& text) {
    std::ofstream out(core::pathFromUtf8(path), std::ios::binary | std::ios::trunc);
    out << text;
    out.close();
    if (!out) {
        throw fs::filesystem_error("cannot write the file", core::pathFromUtf8(path),
                                   std::make_error_code(std::errc::io_error));
    }
}

// write_json: json.dump(indent=2) and a newline, <path>.tmp then os.replace
void writeJson(const std::string& path, const Value& value) {
    const std::string tmp = path + ".tmp";
    writeTextFile(tmp, core::toJsonAscii(value, 2) + "\n");
    fs::rename(core::pathFromUtf8(tmp), core::pathFromUtf8(path));
}

std::string listRepr(const std::vector<std::string>& items) {
    std::string out = "[";
    for (std::size_t k = 0; k < items.size(); ++k) {
        out += (k ? ", " : "") + core::pythonStringRepr(items[k]);
    }
    return out + "]";
}

struct Layout {
    std::string outDir, marker;
    CleanResult result;
};

Layout layoutFor(const core::CleaningParams& params, const std::string& topology,
                 const std::vector<std::string>& trajectories, const std::string& cacheDir) {
    Layout l;
    const std::string key = cleaningKey(params, topology, trajectories);
    l.outDir = abspath(params.output_dir.empty() ? joinPath(joinPath(cacheDir, "cleaned"), key) : params.output_dir);
    CleanResult& r = l.result;
    r.topology = joinPath(l.outDir, isParm7(topology) ? kCleanedParm7 : kCleanedReference);
    r.trajectory = joinPath(l.outDir, kCleanedTrajectory);
    r.referencePdb = joinPath(l.outDir, kCleanedReference);
    r.script = joinPath(l.outDir, "cpptraj.in");
    r.log = joinPath(l.outDir, "cpptraj.log");
    r.outDir = l.outDir;
    r.key = key;
    l.marker = joinPath(l.outDir, kCleaningMarker);
    return l;
}

// the marker's "key" equals `key` (an unreadable marker counts as none)
bool markerMatches(const std::string& marker, const std::string& key) {
    try {
        std::ifstream in(core::pathFromUtf8(marker), std::ios::binary);
        if (!in) {
            return false;
        }
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        const Value info = core::parseJson(text, marker);
        const Value* stored = info.isTable() ? info.find("key") : nullptr;
        return stored && stored->isString() && *stored->get<std::string>() == key;
    } catch (const std::exception&) {
        return false;
    }
}

}  // namespace

std::string findCpptraj(std::string_view pathIn) {
    const std::string path = expanduser(pathIn.empty() ? std::string("cpptraj") : std::string(pathIn));
    if (hasSeparator(path)) {
        if (core::process::isExecutableFile(path)) {
            return abspath(path);
        }
        throw CpptrajError("cpptraj not executable at " + path);
    }
    if (std::optional<std::string> found = core::process::which(path)) {
        return *found;
    }
    const std::string amberhome = core::process::environmentVariable("AMBERHOME");
    if (!amberhome.empty()) {
        const std::string candidate = joinPath(joinPath(amberhome, "bin"), path);
#ifdef _WIN32
        if (exists(candidate)) {
#else
        if (::access(candidate.c_str(), X_OK) == 0) {
#endif
            return candidate;
        }
    }
    throw CpptrajError("cpptraj (" + core::pythonStringRepr(path) +
                       ") was not found on PATH. Install AmberTools or set the cpptraj path in the cleaning options.");
}

std::string buildCpptrajScript(const core::CleaningParams& params, const std::string& topology,
                               const std::vector<std::string>& trajectories, const std::string& outDir) {
    if (trajectories.empty()) {
        throw CpptrajError("cleaning requires at least one trajectory");
    }
    // one push_back per line, never a braced list: quote() throws on a path cpptraj cannot read, and
    // GCC before 13 does not destroy the list's elements made before the throw (GCC bug 66139: the
    // sanitizer job found the leak)
    std::vector<std::string> lines;
    lines.push_back(std::string("# Generated by ") + kName + " " + kVersion + " on " + isoNow());  // CMD-08
    lines.push_back("parm " + quote(abspath(topology)));
    const std::string stop = params.stop && *params.stop != 0 ? std::to_string(*params.stop) : "last";
    for (const std::string& t : trajectories) {
        lines.push_back("trajin " + quote(abspath(t)) + " " + std::to_string(std::max<std::int64_t>(1, params.start)) +
                        " " + stop + " " + std::to_string(std::max<std::int64_t>(1, params.offset)));
    }
    if (params.fit && params.fit_reference == "frame") {
        lines.push_back("reference " + quote(abspath(trajectories.front())) + " " +
                        std::to_string(std::max<std::int64_t>(1, params.reference_frame)) + " [mdatkref]");
    } else if (params.fit && params.fit_reference == "file") {
        if (params.reference_file.empty()) {
            throw CpptrajError("fit reference 'file' selected but no reference file given");
        }
        lines.push_back("reference " + quote(abspath(params.reference_file)) + " [mdatkref]");
    }
    const bool parm7 = isParm7(topology);
    const std::string mask = stripMask(params);
    if (!mask.empty()) {
        lines.push_back("strip " + mask +
                        (parm7 ? " parmout " + quote(joinPath(outDir, kCleanedParm7)) : std::string()));
    } else if (parm7) {
        lines.push_back("parmwrite out " + quote(joinPath(outDir, kCleanedParm7)));
    }
    if (params.autoimage) {
        const std::string anchor = compactMask(params.autoimage_anchor);
        lines.push_back(anchor.empty() ? "autoimage" : "autoimage anchor " + anchor);
    }
    if (params.fit) {
        if (params.remove_box_after_fit) {
            // fitting rotates the unit cell, which NetCDF cannot store
            lines.push_back("box nobox");
        }
        std::string fitMask = compactMask(params.fit_mask);
        if (fitMask.empty()) {
            fitMask = "!@H=";
        }
        const std::string ref = params.fit_reference == "first" ? "first" : "ref [mdatkref]";
        lines.push_back("rms MDATKfit " + fitMask + " " + ref + " mass");
    }
    if (!core::pythonStrip(params.extra_commands).empty()) {
        for (const std::string& l : core::pythonSplitlines(core::pythonStrip(params.extra_commands))) {
            lines.push_back(core::pythonRstrip(l));
        }
    }
    lines.push_back("trajout " + quote(joinPath(outDir, kCleanedTrajectory)) + " netcdf");
    lines.push_back("trajout " + quote(joinPath(outDir, kCleanedReference)) + " pdb onlyframes 1");
    lines.push_back("run");
    lines.push_back("quit");
    std::string script;
    for (const std::string& l : lines) {
        script += l + "\n";
    }
    return script;
}

std::vector<std::string> cpptrajCommand(const core::CleaningParams& params,
                                        const std::function<void(const std::string&)>& log) {
    const auto say = [&log](const std::string& message) {
        if (log) {
            log(message);
        }
    };
    const std::string envSpec = core::pythonStrip(params.conda_env);
    if (!envSpec.empty()) {
        const std::optional<core::conda::CondaEnv> env = core::conda::resolveEnv(envSpec);
        if (const std::optional<std::string> exe = core::conda::findCondaExe(params.conda_exe, env)) {
            return core::conda::runCommand(*exe, env ? env->prefix : envSpec, {"cpptraj"});
        }
        if (env && env->hasCpptraj()) {
            say("No conda executable found; calling " + env->cpptraj() + " directly");
            return {env->cpptraj()};
        }
        throw CpptrajError("conda environment " + core::pythonStringRepr(envSpec) +
                           " not found or has no cpptraj, and no conda/mamba executable is available to activate it");
    }
    try {
        return {findCpptraj(params.cpptraj_path)};
    } catch (const CpptrajError&) {
        if (!params.cpptraj_path.empty() && params.cpptraj_path != "cpptraj") {
            throw;
        }
        const std::optional<core::conda::CondaEnv> env = core::conda::detectAmbertoolsEnv();
        if (!env) {
            throw;
        }
        const std::optional<std::string> exe = core::conda::findCondaExe(params.conda_exe, env);
        say("cpptraj is not on PATH; using conda environment " + core::pythonStringRepr(env->name) + " (" +
            env->prefix + ")");
        if (exe) {
            return core::conda::runCommand(*exe, env->prefix, {"cpptraj"});
        }
        return {env->cpptraj()};
    }
}

CleanResult clean(const core::CleaningParams& params, const std::string& topology,
                  const std::vector<std::string>& trajectories, const std::string& cacheDir,
                  const core::Reporter& reporter, bool force) {
    if (topology.empty() || !exists(topology)) {
        throw CpptrajError("topology not found: " + topology);
    }
    for (const std::string& t : trajectories) {
        if (!exists(t)) {
            throw CpptrajError("trajectory not found: " + t);
        }
    }
    Layout l = layoutFor(params, topology, trajectories, cacheDir);
    CleanResult& result = l.result;
    if (!force && exists(l.marker) && exists(result.topology) && exists(result.trajectory) &&
        markerMatches(l.marker, result.key)) {
        result.reused = true;
        reporter.rawLog("Reusing cleaned trajectory from " + l.outDir);
        return result;
    }

    const std::vector<std::string> command =
        cpptrajCommand(params, [&reporter](const std::string& message) { reporter.rawLog(message); });
    fs::create_directories(core::pathFromUtf8(l.outDir));
    std::error_code ec;
    fs::remove(core::pathFromUtf8(l.marker), ec);
    const std::string script = buildCpptrajScript(params, topology, trajectories, l.outDir);
    writeTextFile(result.script, script);
    reporter.rawLog("cpptraj script written to " + result.script + ":");
    for (const std::string& line : core::pythonSplitlines(script)) {
        reporter.rawLog("    " + line);
    }

    // run_cpptraj: the output into cpptraj.log and, line by line, into the log; percentages as progress
    std::vector<std::string> argv = command;
    argv.push_back("-i");
    argv.push_back(result.script);
    std::ofstream logFile(core::pathFromUtf8(result.log), std::ios::binary | std::ios::trunc);
    if (!logFile) {
        throw fs::filesystem_error("cannot write the file", core::pathFromUtf8(result.log),
                                   std::make_error_code(std::errc::permission_denied));
    }
    core::process::LineSplitter lines;
    core::process::RunOptions options;
    options.workingDirectory = core::pathFromUtf8(l.outDir);
    options.cancelled = [&reporter]() { return reporter.cancelled(); };
    options.onOutput = [&](std::string_view chunk) {
        logFile.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
        logFile.flush();
        if (const std::optional<int> last = lastPercent(chunk)) {
            reporter.progress(static_cast<double>(std::min(100, *last)) / 100.0, "cpptraj");
        }
        for (const std::string& line : lines.add(chunk)) {
            reporter.rawLog(line);
        }
    };
    core::process::Result run;
    try {
        run = core::process::run(argv, options);
    } catch (const core::process::StartError& e) {
        throw CpptrajError(e.what());
    }
    if (const std::optional<std::string> rest = lines.finish()) {
        reporter.rawLog(*rest);
    }
    logFile.close();
    if (run.cancelled) {
        throw core::Cancelled("cpptraj cancelled");
    }
    if (run.exitCode != 0) {
        throw CpptrajError("cpptraj failed with exit code " + std::to_string(run.exitCode) + "; see " + result.log);
    }
    std::vector<std::string> missing;
    for (const std::string* p : {&result.topology, &result.trajectory}) {
        if (!exists(*p)) {
            missing.push_back(*p);
        }
    }
    if (!missing.empty()) {
        throw CpptrajError("cpptraj finished but did not write " + listRepr(missing) + "; see " + result.log);
    }
    core::Array trajectoriesIn;
    for (const std::string& t : trajectories) {
        trajectoriesIn.emplace_back(abspath(t));
    }
    writeJson(l.marker, Table{{"project", "MDATK"},
                              {"mdatk_version", kMdatkReference},
                              {"key", result.key},
                              {"created", isoNow()},
                              {"cpptraj", core::process::commandLineText(command)},
                              {"topology_in", abspath(topology)},
                              {"trajectories_in", std::move(trajectoriesIn)},
                              {"params", paramsTable(params)},
                              {"outputs", Table{{"topology", result.topology},
                                                {"trajectory", result.trajectory},
                                                {"reference_pdb", result.referencePdb}}}});
    reporter.progress(1.0, "cpptraj");
    return result;
}

std::optional<CleanResult> cleanedIfAvailable(const core::CleaningParams& params, const std::string& topology,
                                              const std::vector<std::string>& trajectories,
                                              const std::string& cacheDir) {
    try {
        Layout l = layoutFor(params, topology, trajectories, cacheDir);
        if (!exists(l.result.topology)) {  // MDATK: the parm7, else the reference PDB
            l.result.topology = l.result.referencePdb;
        }
        if (exists(l.marker) && exists(l.result.topology) && exists(l.result.trajectory) &&
            markerMatches(l.marker, l.result.key)) {
            l.result.reused = true;
            return l.result;
        }
    } catch (const std::exception&) {  // OSError: an input is missing
    }
    return std::nullopt;
}

}  // namespace mdbridge::analysis
