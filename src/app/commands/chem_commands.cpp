#include "chem_commands.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "command_support.hpp"
#include "density_pipeline.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/chem/annotate.hpp"
#include "mdbridge/chem/density_fit.hpp"
#include "mdbridge/chem/dictionaries.hpp"
#include "mdbridge/chem/fit_to_map.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/mutate.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/cif.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/molecule/density_shadow.hpp"
#include "render_pipeline.hpp"

namespace mdbridge::app::detail {

namespace fs = std::filesystem;
using core::Value;

// ---- the parsed options -----------------------------------------------------------------------

struct ChemCommands::Arguments {
    const NetworkBackend* network = nullptr;

    // pdb
    CLI::App* pdb = nullptr;
    CLI::App* clean = nullptr;
    CLI::App* check = nullptr;
    CLI::App* toLrna = nullptr;
    CLI::App* annotate = nullptr;
    std::string pdbProject;
    std::string input;
    std::string output;
    std::string scheme, records;
    bool noTerminalNames = false, keepConect = false, noDetect = false;
    std::string lRanges, chainReference, chainIds, splitAt, combine, combineInto;
    CLI::Option* combineIntoOption = nullptr;
    bool assignChains = false, copyNumbers = false, splitSegments = false, sortChains = false;
    bool renumber = false, perChain = false, keepGaps = false, seqres = false, strip5p = false;
    std::int64_t start = 1;
    std::string fasta;
    std::string ranges;
    std::string toLrnaScheme = "amber";
    bool noPairs = false, noStacking = false, noDotbracket = false;
    std::string engine = "auto";

    // build
    CLI::App* build = nullptr;
    CLI::App* helix = nullptr;
    CLI::App* mutate = nullptr;
    CLI::App* fit = nullptr;
    std::string buildProject;
    std::string seq, partner = "auto", chirality = "D", pairing = "wc", triplex = "none", third = "auto";
    std::string chains = "A,B", startNumbers = "1,1", helixScheme;
    bool helixNoTerminalNames = false, phosphate5 = false, hydrogens = false;
    std::string buildOutput;
    std::string buildInput, mutations, map;
    bool noCompensatory = false, wobble = false;
    std::optional<double> fitResolution;
    std::string fitSearch = "global";
    double fitStep = 30.0;

    // dict
    CLI::App* dict = nullptr;
    CLI::App* dictList = nullptr;
    CLI::App* dictSmiles = nullptr;
    CLI::App* dictFinish = nullptr;
    CLI::App* dictPubchem = nullptr;
    bool listJson = false;
    std::string code;
    std::vector<std::string> finishInputs;
    std::string finishOutput, group = std::string(chem::kDefaultGroup), name, sourceId, finishProject;
    bool allowProtonation = false, keepNames = false;
    std::string compound, pubchemCode;
    double timeout = 15.0;

    // density fit / simulate
    CLI::App* densityFit = nullptr;
    CLI::App* densitySimulate = nullptr;
    CommonArguments fitCommon;
    std::string fitMap;
    std::optional<double> resolution, level, step, spacing;
    std::string frame;
    CLI::Option* frameOption = nullptr;
    std::string search;
    CLI::Option* saveOption = nullptr;
    std::string save;
    CommonArguments simulateCommon;
    std::string simulateFrame;
    bool placed = false;
    std::string simulateOutput;
};

namespace {

// ---- options with Python's readings ----------------------------------------------------------

// type=float: float(text) (also nan, inf, 1_000, surrounding spaces); argparse's usage error else
CLI::Option* floatOption(CLI::App& command, const std::string& names, std::optional<double>& target,
                         const std::string& description) {
    CLI::Option* option = command.add_option_function<std::string>(
        names, [&target](const std::string& text) { target = core::pythonFloat(text); }, description);
    option->check(CLI::Validator(
        [](std::string& text) { return core::tryPythonFloat(text) ? std::string() : "invalid float value: " + text; },
        "FLOAT"));
    return lastWins(option);
}

CLI::Option* floatOption(CLI::App& command, const std::string& names, double& target, const std::string& description) {
    CLI::Option* option = command.add_option_function<std::string>(
        names, [&target](const std::string& text) { target = core::pythonFloat(text); }, description);
    option->check(CLI::Validator(
        [](std::string& text) { return core::tryPythonFloat(text) ? std::string() : "invalid float value: " + text; },
        "FLOAT"));
    return lastWins(option);
}

// type=int: int(text)
CLI::Option* intOption(CLI::App& command, const std::string& names, std::int64_t& target,
                       const std::string& description) {
    CLI::Option* option = command.add_option_function<std::string>(
        names, [&target](const std::string& text) { target = core::pythonInt(text); }, description);
    option->check(CLI::Validator(
        [](std::string& text) { return core::tryPythonInt(text) ? std::string() : "invalid int value: " + text; },
        "INT"));
    return lastWins(option);
}

CLI::Option* choice(CLI::Option* option, std::vector<std::string> choices) {
    return option->check(CLI::IsMember(std::move(choices)));
}

std::vector<std::string> toStrings(std::span<const std::string_view> items) {
    return {items.begin(), items.end()};
}

// ---- files ------------------------------------------------------------------------------------

// open(path, "w", encoding="utf-8").write(text): OSError texts name the path as given
void writeText(const std::string& path, const std::string& text) {
    std::ofstream out(core::pathFromUtf8(path), std::ios::binary | std::ios::trunc);
    if (!out) {
        const int error = errno != 0 ? errno : EIO;
        std::error_code code(error, std::generic_category());
        std::error_code ec;
        if (fs::is_directory(core::pathFromUtf8(path), ec)) {
            code = std::make_error_code(std::errc::is_a_directory);
        }
        throw fs::filesystem_error("open", core::pathFromUtf8(path), code);
    }
    out << text;
    out.close();
    if (!out) {
        throw fs::filesystem_error("write", core::pathFromUtf8(path), std::make_error_code(std::errc::io_error));
    }
}

// open(path, encoding="utf-8", errors="replace").read(): bytes that are not UTF-8 as U+FFFD,
// universal newlines
std::string readTextReplace(const std::string& path) {
    std::ifstream in(core::pathFromUtf8(path), std::ios::binary);
    if (!in) {
        std::error_code ec;
        const std::error_code code = fs::is_directory(core::pathFromUtf8(path), ec)
                                         ? std::make_error_code(std::errc::is_a_directory)
                                         : std::error_code(errno != 0 ? errno : ENOENT, std::generic_category());
        throw fs::filesystem_error("open", core::pathFromUtf8(path), code);
    }
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string text = core::encodeUtf8(core::decodeUtf8(bytes));
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '\r') {
            out += '\n';
            if (i + 1 < text.size() && text[i + 1] == '\n') ++i;
        } else {
            out += text[i];
        }
    }
    return out;
}

// CH-13: a PDB input is a file; a missing one is named ([Errno 2] ...), a folder too ([Errno 21])
chem::PdbSource pdbFile(const std::string& path) {
    std::error_code ec;
    const fs::path p = core::pathFromUtf8(path);
    if (fs::is_directory(p, ec)) {
        throw fs::filesystem_error("open", p, std::make_error_code(std::errc::is_a_directory));
    }
    if (!fs::exists(p, ec)) {
        throw fs::filesystem_error("open", p, std::make_error_code(std::errc::no_such_file_or_directory));
    }
    return chem::PdbSource::file(p);
}

std::string joined(const std::vector<std::string>& lines, const char* separator = "\n") {
    std::string out;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        out += (i ? separator : "") + lines[i];
    }
    return out;
}

void printNotes(Console& console, const std::string& prefix, const std::vector<std::string>& notes) {
    for (const std::string& note : notes) {
        console.err << prefix << ": " << note << '\n';
    }
    console.err.flush();
}

// [c.strip() for c in text.split(",") if c.strip()] if text else None
std::optional<std::vector<std::string>> ids(const std::string& text) {
    if (text.empty()) {
        return std::nullopt;
    }
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = text.find(',', start);
        const std::string part = core::pythonStrip(
            std::string_view(text).substr(start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!part.empty()) out.push_back(part);
        if (comma == std::string::npos) return out;
        start = comma + 1;
    }
}

chem::RangeSpec rangesOrNone(const std::string& text) {
    return text.empty() ? chem::RangeSpec() : chem::RangeSpec(text);
}

std::vector<std::string> keysFor(std::initializer_list<std::string> resources) {
    const std::vector<std::string> ids(resources);
    return core::citations::keysFor(ids);
}

// _project_file(args): the project named by --project, as given (core::logged skips parameter logs)
fs::path projectFile(const std::string& params) { return params.empty() ? fs::path() : core::pathFromUtf8(params); }

// ---- pdb ----------------------------------------------------------------------------------------

int pdbRun(ChemCommands::Arguments& a, Console& console, const std::vector<std::string>& arguments) {
    if (a.check->parsed()) {  // pdbtools.run: no step log
        const Value info = chem::check(pdbFile(a.input), rangesOrNone(a.lRanges));
        console.out << core::toJsonAscii(info, 2) << '\n';
        return 0;
    }
    const std::string name = core::pathBasename(a.input);
    if (a.annotate->parsed()) {
        core::logged(projectFile(a.pdbProject), "pdb annotate", "Base-pair annotation of " + name, {},
                     [&](core::LogEntry& entry) {
                         entry.command = commandText(arguments, "pdb");
                         const chem::PdbSource source = pdbFile(a.input);  // CH-13
                         chem::AnnotateOptions options;
                         options.pairs = !a.noPairs;
                         options.stacking = !a.noStacking;
                         options.dotbracket = !a.noDotbracket;
                         options.engine = a.engine;
                         const chem::AnnotationResult res = chem::annotate(source, options);
                         writeText(a.output, res.text);
                         printNotes(console, "annotate", res.notes);
                         entry.files = {a.output};
                         entry.summary = joined(res.notes);
                         // CH-14: never Barnaba, so never its citation
                     });
        return 0;
    }
    const std::string command = a.toLrna->parsed() ? "to-lrna" : "clean";
    core::logged(projectFile(a.pdbProject), "pdb " + command, "PDB clean-up of " + name, {},
                 [&](core::LogEntry& entry) {
                     entry.command = commandText(arguments, "pdb");
                     chem::CleanResult res;
                     if (a.toLrna->parsed()) {
                         res = chem::toLrna(pdbFile(a.input), chem::RangeSpec(a.ranges), a.toLrnaScheme, true,
                                            !a.noTerminalNames);
                     } else {
                         chem::CleanOptions o;
                         if (!a.scheme.empty()) o.scheme = a.scheme;
                         if (!a.records.empty()) o.records = a.records;
                         o.stripConect = !a.keepConect;
                         o.detectL = !a.noDetect;
                         o.lRanges = rangesOrNone(a.lRanges);
                         const chem::PdbSource source = pdbFile(a.input);
                         if (!a.chainReference.empty()) o.chainReference = pdbFile(a.chainReference);
                         o.chainIds = ids(a.chainIds);
                         o.assignChains = a.assignChains;
                         o.copyNumbers = a.copyNumbers;
                         o.splitSegments = a.splitSegments;
                         o.splitAt = rangesOrNone(a.splitAt);
                         o.combine = ids(a.combine).value_or(std::vector<std::string>{});
                         if (a.combineIntoOption->count() > 0) o.combineInto = a.combineInto;
                         o.sortChains = a.sortChains;
                         o.renumber = a.renumber;
                         o.start = a.start;
                         o.perChain = a.perChain;
                         o.keepGaps = a.keepGaps;
                         o.seqres = a.seqres;
                         o.makeFasta = !a.fasta.empty();
                         o.terminalNames = !a.noTerminalNames;
                         o.strip5p = a.strip5p;
                         res = chem::clean(source, o);
                     }
                     if (!a.output.empty()) {
                         writeText(a.output, res.text);
                     } else {
                         console.out << res.text;
                         console.out.flush();
                     }
                     if (a.clean->parsed() && !a.fasta.empty() && res.fasta && !res.fasta->empty()) {
                         writeText(a.fasta, *res.fasta);
                     }
                     printNotes(console, "pdbtools", res.notes);
                     entry.files.clear();
                     if (!a.output.empty()) entry.files.push_back(a.output);
                     if (a.clean->parsed() && !a.fasta.empty()) entry.files.push_back(a.fasta);
                     entry.summary = joined(res.notes);
                 });
    return 0;
}

// ---- build --------------------------------------------------------------------------------------

int buildRun(ChemCommands::Arguments& a, Console& console, const std::vector<std::string>& arguments) {
    std::string step;
    std::string title;
    if (a.helix->parsed()) {
        step = "build helix";
        title = "Ideal " + chem::helixKind(a.pairing, a.triplex) + " " + a.chirality + "-RNA helix " +
                core::pythonStrip(a.seq);
    } else if (a.mutate->parsed()) {
        step = "build mutate";
        title = "Mutation of " + core::pathBasename(a.buildInput) + ": " + core::pythonStrip(a.mutations);
    } else {
        step = "build fit";
        title = "Rigid-body fit of " + core::pathBasename(a.buildInput) + " into " + core::pathBasename(a.map);
    }
    std::optional<double> cc;
    core::logged(projectFile(a.buildProject), step, title, {}, [&](core::LogEntry& entry) {
        entry.command = commandText(arguments, "build");
        std::vector<std::string> notes;
        std::string text;
        if (a.helix->parsed()) {
            chem::HelixOptions o;
            o.partner = a.partner;
            o.chirality = a.chirality;
            o.chains = a.chains;
            o.start = a.startNumbers;
            if (!a.helixScheme.empty()) o.scheme = a.helixScheme;
            o.terminals = !a.helixNoTerminalNames;
            o.hydrogens = a.hydrogens;
            o.phosphate5 = a.phosphate5;
            o.pairing = a.pairing;
            o.triplex = a.triplex;
            o.third = a.third;
            text = chem::buildHelix(a.seq, o, &notes);
        } else if (a.mutate->parsed()) {
            chem::MutateOptions o;
            o.compensatory = !a.noCompensatory;
            o.wobble = a.wobble;
            chem::MutateResult res = chem::mutate(pdbFile(a.buildInput), a.mutations, o);  // CH-13
            text = std::move(res.text);
            notes = std::move(res.notes);
        } else {
            chem::FitToMapOptions o;
            o.search = a.fitSearch;
            o.angularStep = a.fitStep;
            chem::FitToMapResult res = chem::fitToMap(pdbFile(a.buildInput), core::pathFromUtf8(a.map),
                                                      a.fitResolution.value_or(0.0), o);
            text = std::move(res.text);
            notes = std::move(res.notes);
            cc = res.cc;
        }
        writeText(a.buildOutput, text);
        printNotes(console, "build", notes);
        entry.files = {a.buildOutput};
        entry.summary = joined(notes);
        if (a.helix->parsed()) {  // the template sources
            const bool other = a.pairing != "wc" || !(a.triplex.empty() || a.triplex == "none");
            // D19: the Hoogsteen duplex and the pyrimidine-motif triplex also cite their RNA triple helices
            const bool triplexes = a.pairing == "hoogsteen" || (a.pairing == "wc" && a.triplex == "hoogsteen");
            entry.citations = triplexes ? keysFor({"arnott-arna", "pdb-4jrd", "pdb-triplexes"})
                              : other   ? keysFor({"arnott-arna", "pdb-4jrd"})
                                        : keysFor({"arnott-arna"});
        }
    });
    if (cc) {
        console.out << core::toJsonAscii(Value(core::Table{{"output", absolutePath(a.buildOutput)},
                                                           {"cc", core::pythonRound(*cc, 4)}}),
                                         2)
                    << '\n';
    }
    return 0;
}

// ---- dict ---------------------------------------------------------------------------------------

int dictRun(ChemCommands::Arguments& a, Console& console, const std::vector<std::string>& arguments) {
    if (a.dictList->parsed()) {
        core::Array rows;
        for (const chem::Component& c : chem::dictionaryComponents()) {
            core::Array leaving;
            for (const std::string& l : c.leaving) leaving.emplace_back(l);
            rows.emplace_back(core::Table{{"code", c.code},
                                          {"parent", c.parent},
                                          {"terminal", c.terminal},
                                          {"name", c.name},
                                          {"smiles", c.smiles},
                                          {"leaving", std::move(leaving)},
                                          {"formula", c.formula()}});
        }
        if (a.listJson) {
            console.out << core::toJsonAscii(Value(std::move(rows)), 2) << '\n';
        } else {
            for (const chem::Component& c : chem::dictionaryComponents()) {
                std::string code = c.code;
                while (core::characterCount(code) < 4) code += ' ';
                console.out << code << ' ' << c.smiles
                            << (c.leaving.empty() ? std::string()
                                                  : "  (" + joined(c.leaving, ", ") + " removed after generation)")
                            << '\n';
            }
        }
        return 0;
    }
    if (a.dictSmiles->parsed()) {
        const chem::Component& comp = chem::getComponent(a.code);
        console.out << comp.smiles << '\n';
        console.out.flush();
        if (!comp.leaving.empty()) {
            console.err << "dictionaries: " << comp.code << " is generated from this molecule; "
                        << joined(comp.leaving, ", ") << " is removed by 'finish'\n";
        }
        return 0;
    }
    if (a.dictFinish->parsed()) {
        core::logged(projectFile(a.finishProject), "dict finish", "Restraint dictionary " + a.code, keysFor({"ccd"}),
                     [&](core::LogEntry& entry) {
                         entry.command = commandText(arguments, "dict");
                         const std::vector<std::string> exclude{absolutePath(a.finishOutput)};
                         const std::string src = io::findDictionary(a.finishInputs, a.code, exclude);
                         const std::string text = readTextReplace(src);
                         chem::FinishOptions o;
                         o.group = a.group;
                         o.name = a.name;
                         o.sourceId = a.sourceId;
                         o.allowProtonation = a.allowProtonation;
                         o.keepNames = a.keepNames;
                         const chem::FinishResult res = chem::finishDictionary(text, a.code, o);
                         writeText(a.finishOutput, res.text);
                         std::vector<std::string> notes{src + " (" +
                                                        (res.program.empty() ? std::string("unknown program") : res.program) +
                                                        ", component " + res.sourceId + ") -> " + a.finishOutput};
                         notes.insert(notes.end(), res.notes.begin(), res.notes.end());
                         printNotes(console, "dictionaries", notes);
                         entry.files = {a.finishOutput};
                         entry.summary = joined(notes);
                     });
        return 0;
    }
    // pubchem (network, on request only: owner decision D7)
    const std::string smiles =
        fetchPubchemSmiles(a.compound, a.timeout, a.network ? a.network->get : HttpGet());
    console.out << smiles << '\n';
    console.out.flush();
    std::vector<std::string> codes;
    if (!a.pubchemCode.empty()) {
        std::string upper;
        for (const char32_t c : core::decodeUtf8Surrogateescape(a.pubchemCode)) {
            core::appendUtf8(upper, (c >= U'a' && c <= U'z') ? c - 32 : c);
        }
        codes.push_back(upper);
    } else {
        for (const chem::Component& c : chem::dictionaryComponents()) codes.push_back(c.code);
    }
    std::vector<std::string> hits;
    std::vector<std::string> problems;
    for (const std::string& code : codes) {
        std::vector<std::string> notes;
        try {
            notes = chem::checkSmiles(code, smiles);
        } catch (const io::DictionaryError& error) {
            problems.push_back(code + ": " + error.what());
            continue;
        }
        hits.push_back(code);
        for (const std::string& n : notes) {
            console.err << "dictionaries: " << code << ": " << n << '\n';
        }
    }
    if (!hits.empty()) {
        console.err << "dictionaries: matches " << joined(hits, ", ") << '\n';
    } else {
        console.err << "dictionaries: matches no known code"
                    << (!a.pubchemCode.empty() && !problems.empty() ? " (" + problems.front() + ")" : std::string())
                    << '\n';
    }
    return 0;
}

// ---- density fit / simulate ---------------------------------------------------------------------

bool positiveFinite(double x) { return std::isfinite(x) && x > 0; }

// CH-12: the project file as written (ProjectDocument: paths as the file has them), with the fit's
// density settings; a --map given relative stored relative to the project's folder; --save OTHER
// rebases the relative paths onto the other folder (Save As)
void saveAsWritten(const ChemCommands::Arguments& a, const core::DensityParams& fitted, const std::string& target) {
    const std::string source = absolutePath(a.fitCommon.params);
    core::ProjectDocument document;
    document.open(core::pathFromUtf8(source));
    document.edit("Density fit", [&](core::Project& written) {
        core::DensityParams& w = written.density;
        if (!a.fitMap.empty()) {
            if (core::pathFromUtf8(a.fitMap).is_absolute()) {
                w.map_path = a.fitMap;
            } else {
                const std::optional<std::string> rel = core::relativePath(
                    absolutePath(a.fitMap), core::pathToUtf8(core::pathFromUtf8(source).parent_path()));
                w.map_path = rel ? *rel : absolutePath(a.fitMap);
            }
        }
        w.resolution = fitted.resolution;
        w.level = fitted.level;
        w.fit_frame = fitted.fit_frame;
        w.angular_step = fitted.angular_step;
        w.transform = fitted.transform;
        w.anchor = fitted.anchor;
        w.fit_cc = fitted.fit_cc;
        w.fit_info = fitted.fit_info;
    });
    const std::string targetAbs = absolutePath(target);
    if (targetAbs == source && !document.isParamsLog()) {
        document.save();
    } else {
        document.saveAs(core::pathFromUtf8(target));
    }
}

int densityFitRun(ChemCommands::Arguments& a, CLI::App& command, Console& console) {
    if (a.resolution && !std::isfinite(*a.resolution)) {  // CH-10
        throw core::ValueError("the resolution must be a positive number of Å");
    }
    core::Project project = loadProject(command, a.fitCommon, true);
    // _density_overrides
    core::DensityParams& d = project.density;
    if (!a.fitMap.empty()) {
        const std::string map = absolutePath(a.fitMap);
        if (map != d.map_path) {
            d.transform.clear();
            d.fit_cc = 0.0;
            d.fit_info.clear();
        }
        d.map_path = map;
    }
    if (a.resolution) d.resolution = *a.resolution;
    if (a.level) d.level = *a.level;
    if (a.frameOption->count() > 0) d.fit_frame = a.frame;
    if (a.step) d.angular_step = *a.step;
    std::vector<std::string> problems = core::densityProblems(project);
    if (d.map_path.empty()) {
        problems.push_back("no density map: give --map or set density.map_path in the project");
    }
    if (!problems.empty()) {
        throw core::ValueError(joined(problems, "; "));
    }
    ConsoleReporter consoleReporter(console, a.fitCommon.quiet, a.fitCommon.verbose);
    analysis::Pipeline pipeline(std::move(project), consoleReporter.reporter());
    checkProject(pipeline.project(), pipeline.reporter());
    const chem::FitResult res =
        fitDensity(pipeline, a.search.empty() ? std::nullopt : std::optional<std::string>(a.search));
    const core::DensityParams& fitted = pipeline.project().density;
    core::Array transform;
    for (const double v : fitted.transform) transform.emplace_back(v);
    core::Table out{{"cc", core::pythonRound(res.cc, 4)},
                    {"start", res.start},
                    {"info", fitted.fit_info},
                    {"transform", std::move(transform)}};
    if (a.saveOption->count() > 0) {
        const std::string target = !a.save.empty() ? a.save : a.fitCommon.params;
        if (target.empty()) {
            throw core::ValueError("--save needs a path when no --project was given");
        }
        if (!a.fitCommon.params.empty()) {
            saveAsWritten(a, fitted, target);  // CH-12
        } else {
            core::saveProject(core::pathFromUtf8(target), pipeline.project());
        }
        out.emplace_back("saved", absolutePath(target));
    }
    console.out << core::toJsonAscii(Value(std::move(out)), 2) << '\n';
    return 0;
}

int densitySimulateRun(ChemCommands::Arguments& a, CLI::App& command, Console& console) {
    const double resolution = a.resolution.value_or(0.0);
    if (!positiveFinite(resolution)) {  // CH-10
        throw core::ValueError("the resolution must be a positive number of Å");
    }
    if (a.spacing && !(std::isfinite(*a.spacing) && *a.spacing >= 0)) {
        throw core::ValueError("the grid spacing must be positive");
    }
    core::Project project = loadProject(command, a.simulateCommon, true);
    ConsoleReporter consoleReporter(console, a.simulateCommon.quiet, a.simulateCommon.verbose);
    analysis::Pipeline pipeline(std::move(project), consoleReporter.reporter());
    checkProject(pipeline.project(), pipeline.reporter());
    const DensityFrameAtoms atoms = densityFrameAtoms(pipeline, a.simulateFrame.empty() ? "ref" : a.simulateFrame);
    std::vector<core::Vec3> xyz = atoms.xyz;
    const core::DensityParams& d = pipeline.project().density;
    if (a.placed) {
        if (d.transform.size() != 12) {
            throw core::ValueError("--placed needs a stored placement (run 'mdbridge density fit --save' first)");
        }
        xyz = render::molecule::anchored(d.transform, d.anchor, atoms.reference).transform.apply(xyz);
    }
    const std::optional<double> spacing =
        a.spacing && *a.spacing != 0.0 ? a.spacing : std::nullopt;  // `args.spacing or None`
    const io::DensityMap map = io::simulateMap(xyz, io::atomWeights(atoms.elements), resolution, spacing);
    io::writeMap(core::pathFromUtf8(a.simulateOutput), map);
    core::Array grid;
    for (const std::int64_t n : map.shapeXyz) grid.emplace_back(n);
    console.out << core::toJsonAscii(
                       Value(core::Table{{"map", absolutePath(a.simulateOutput)},
                                         {"frame", atoms.frame},
                                         {"atoms", static_cast<std::int64_t>(xyz.size())},
                                         {"grid_xyz", std::move(grid)},
                                         {"voxel_A", core::pythonRound(map.voxel[0], 4)},
                                         {"frame_of_reference", a.placed ? "map (placed)" : "model (superposed)"}}),
                       2)
                << '\n';
    return 0;
}

}  // namespace

// ---- the command line -------------------------------------------------------------------------

ChemCommands::ChemCommands(CLI::App& app, CLI::App& density, const NetworkBackend* network)
    : m_args(std::make_unique<Arguments>()) {
    Arguments& a = *m_args;
    a.network = network;
    const std::vector<std::string> schemes = toStrings(chem::kSchemes);

    // ---- pdb (pdbtools.add_commands, cli.py)
    a.pdb = app.add_subcommand("pdb", "PDB hygiene: ChimeraX repair, L-RNA detection and naming (Amber / CCD / "
                                      "Rosetta), chain ids, renumbering, SEQRES / FASTA");
    a.pdb->require_subcommand(0, 1);
    a.clean = a.pdb->add_subcommand("clean", "repair, rename, re-chain, renumber and rewrite a PDB file");
    {
        CLI::App& c = *a.clean;
        c.add_option("input", a.input, "PDB file")->required();
        lastWins(c.add_option("--output,-o", a.output, "output PDB (default: print to stdout)"));
        choice(lastWins(c.add_option("--scheme", a.scheme,
                                     "residue names: amber (LA/LA5... for the L-RNA force field), ccd (0A... for "
                                     "Coot/Phenix), rosetta (0A..., all ATOM), pdb")),
               schemes);
        choice(lastWins(c.add_option("--records", a.records, "ATOM/HETATM policy (default: follows the scheme)")),
               toStrings(chem::kRecordPolicies));
        c.add_flag("--no-terminal-names", a.noTerminalNames,
                   "chain ends keep the internal names (LA, not LA5/LA3), for libraries that map terminals with "
                   "addPdbResMap");
        c.add_flag("--keep-conect", a.keepConect,
                   "keep CONECT records (stripped by default until final deposition processing)");
        c.add_flag("--no-detect", a.noDetect,
                   "do not detect L-nucleotides from the sugar chirality (C1'/C3'/C4'); use the residue names");
        lastWins(c.add_option("--l-ranges", a.lRanges, "force L for these residues, e.g. 'C:45-73' or 'C' (whole chain)"));
        lastWins(c.add_option("--chain-reference", a.chainReference,
                              "copy chain ids residue by residue from this PDB (e.g. the structure a cpptraj PDB was "
                              "built from)"));
        lastWins(c.add_option("--chain-ids", a.chainIds, "chain ids per TER segment, e.g. A,B,C"));
        c.add_flag("--assign-chains", a.assignChains, "assign A, B, C ... per TER segment");
        c.add_flag("--copy-numbers", a.copyNumbers, "with --chain-reference: also copy the residue numbers");
        c.add_flag("--split-segments", a.splitSegments, "give every TER segment its own new chain id");
        lastWins(c.add_option("--split-at", a.splitAt, "split a chain before these residues, e.g. 'A:31'"));
        lastWins(c.add_option("--combine", a.combine, "merge these chains, e.g. A,B"));
        a.combineIntoOption =
            lastWins(c.add_option("--combine-into", a.combineInto, "chain id of the merged chain (default: the first)"));
        c.add_flag("--sort-chains", a.sortChains, "write the chains in alphabetical order");
        c.add_flag("--renumber", a.renumber, "renumber residues continuously across chains");
        intOption(c, "--start", a.start, "first residue number (default 1)");
        c.add_flag("--per-chain", a.perChain, "restart the numbering in every chain");
        c.add_flag("--keep-gaps", a.keepGaps, "keep numbering gaps when renumbering");
        c.add_flag("--seqres", a.seqres, "write SEQRES records");
        c.add_flag("--strip-5p", a.strip5p,
                   "remove the phosphate (P/OP1/OP2/OP3) of every chain's 5'-terminal nucleotide, as Amber's 5' "
                   "units expect");
        lastWins(c.add_option("--fasta", a.fasta, "also write the sequences to this FASTA file"));
    }
    a.check = a.pdb->add_subcommand("check", "report chains, sequences, residue names and sugar chirality");
    a.check->add_option("input", a.input, "PDB file")->required();
    lastWins(a.check->add_option("--l-ranges", a.lRanges, "residues to treat as L, e.g. 'C:45-73'"));
    a.toLrna = a.pdb->add_subcommand("to-lrna", "rename residues in ranges to L-RNA names (replaces to_lrna.sh)");
    a.toLrna->add_option("input", a.input, "PDB file")->required();
    a.toLrna->add_option("ranges", a.ranges, "e.g. 'C:45-73' or 'A:1-5,A:9-12'")->required();
    choice(lastWins(a.toLrna->add_option("--scheme", a.toLrnaScheme, "residue names (default amber)")), schemes);
    a.toLrna->add_flag("--no-terminal-names", a.noTerminalNames, "plain LG/LA/LC/LU at chain ends, like to_lrna.sh");
    lastWins(a.toLrna->add_option("--output,-o", a.output, "output PDB (default: stdout)"));
    for (CLI::App* c : {a.clean, a.check, a.toLrna}) {
        lastWins(c->add_option("--project,-p,--params", a.pdbProject, "project file whose step log records the clean-up"));
    }
    a.annotate = a.pdb->add_subcommand(
        "annotate", "add REMARK 950 base-pair / dot-bracket records (MDBridge's builtin canonical-pair detector)");
    a.annotate->add_option("input", a.input, "PDB file")->required();
    lastWins(a.annotate->add_option("--output,-o", a.output, "output PDB"))->required();
    a.annotate->add_flag("--no-pairs", a.noPairs, "no PAIR records");
    a.annotate->add_flag("--no-stacking", a.noStacking, "no STACK records");
    a.annotate->add_flag("--no-dotbracket", a.noDotbracket, "no sequence / dot-bracket records");
    choice(lastWins(a.annotate->add_option("--engine", a.engine,
                                           "auto and builtin: the builtin detector; barnaba: the builtin detector "
                                           "with a warning (MDBridge has no Barnaba)")),
           toStrings(chem::kAnnotationEngines));
    lastWins(a.annotate->add_option("--project,-p,--params", a.pdbProject,
                                    "project file whose step log records the annotation"));

    // ---- build (builder.add_commands)
    a.build = app.add_subcommand("build", "starting structures: ideal A-form helices (D / L), base mutations, "
                                          "rigid-body fit into a density map");
    a.build->require_subcommand(0, 1);
    a.helix = a.build->add_subcommand("helix", "ideal A-form RNA helix (D or L), single strand or duplex");
    {
        CLI::App& h = *a.helix;
        lastWins(h.add_option("--seq,--sequence", a.seq, "strand 1, 5'->3' (e.g. GGCAU)"))->required();
        lastWins(h.add_option("--partner", a.partner,
                              "auto (Watson-Crick complement), none (single strand) or the partner sequence 5'->3' "
                              "(same length; mismatches / G·U allowed)"));
        choice(lastWins(h.add_option("--chirality", a.chirality, "D (right-handed) or L (mirror image)")), {"D", "L"});
        choice(lastWins(h.add_option("--pairing", a.pairing,
                                     "wc: antiparallel Watson-Crick A-form (default); hoogsteen (cis WC/Hoogsteen, "
                                     "strand 1 the purine strand: A and G only) and hh (trans Hoogsteen/Hoogsteen "
                                     "A·A, as PDB 4JRD: A only): parallel strands. reverse-wc (trans WC/WC, "
                                     "parallel) and "
                                     "reverse-hoogsteen (trans WC/Hoogsteen, antiparallel) are refused: RNA is not "
                                     "known to form them (DNA and DNA-RNA hybrid builds will come with DNA structure "
                                     "building)")),
               toStrings(chem::kPairings));
        std::vector<std::string> triplexes{"none"};
        for (const std::string_view m : chem::kTriplexMotifs) triplexes.emplace_back(m);
        choice(lastWins(h.add_option("--triplex", a.triplex,
                                     "add a third strand on the Hoogsteen edge of the strand-1 purines (strand 1 A "
                                     "and G only): hoogsteen (pyrimidine motif, parallel); needs --pairing wc. "
                                     "reverse-hoogsteen (purine "
                                     "motif, antiparallel) is refused: RNA is not known to form it (DNA builds will "
                                     "come with DNA structure building)")),
               triplexes);
        lastWins(h.add_option("--third-strand,--third", a.third,
                              "third strand 5'->3' (same length) or auto (the canonical motif bases)"));
        lastWins(h.add_option("--chains", a.chains,
                              "chain ids of strand 1 and 2 (and 3 of a triplex; default A,B / A,B,C)"));
        lastWins(h.add_option("--start", a.startNumbers,
                              "first residue numbers of strand 1 and 2 (and 3; default 1,1; write negative numbers as "
                              "--start=-5,1)"));
        choice(lastWins(h.add_option("--scheme", a.helixScheme,
                                     "residue names (default: A/C/G/U and 0A/0C/0G/0U, no terminal variants)")),
               schemes);
        h.add_flag("--no-terminal-names", a.helixNoTerminalNames, "no 5'/3' residue name variants");
        h.add_flag("--phosphate5", a.phosphate5, "keep a 5'-phosphate (P, OP1, OP2, OP3)");
        h.add_flag("--hydrogens", a.hydrogens, "add all hydrogens (default: heavy atoms only)");
        lastWins(h.add_option("--output,-o", a.buildOutput, "output PDB"))->required();
    }
    a.mutate = a.build->add_subcommand("mutate", "swap bases keeping the glycosidic frame (with compensatory partners)");
    a.mutate->add_option("input", a.buildInput, "PDB file")->required();
    a.mutate->add_option("mutations", a.mutations, "e.g. 'A:5G, B:20C'")->required();
    a.mutate->add_flag("--no-compensatory", a.noCompensatory, "do not change the paired partners");
    a.mutate->add_flag("--wobble", a.wobble, "keep G·U pairs where the new base allows them");
    lastWins(a.mutate->add_option("--output,-o", a.buildOutput, "output PDB"))->required();
    a.fit = a.build->add_subcommand("fit", "rigid-body fit of a model into a density map (MRC / CCP4)");
    a.fit->add_option("input", a.buildInput, "model PDB")->required();
    a.fit->add_option("map", a.map, "density map")->required();
    floatOption(*a.fit, "--resolution", a.fitResolution, "map resolution in Å")->required();
    choice(lastWins(a.fit->add_option("--search", a.fitSearch,
                                      "global search (default) or local refinement of the model as placed")),
           {"global", "local"});
    floatOption(*a.fit, "--step", a.fitStep, "angular step of the global search (degrees)");
    lastWins(a.fit->add_option("--output,-o", a.buildOutput, "output PDB"))->required();
    for (CLI::App* c : {a.helix, a.mutate, a.fit}) {
        lastWins(c->add_option("--project,-p,--params", a.buildProject, "project file whose step log records the step"));
    }

    // ---- dict (dictionaries.add_commands)
    a.dict = app.add_subcommand("dict", "restraint dictionaries for L-RNA residues: SMILES for acedrg / eLBOW, "
                                        "renaming to the CCD atom names, PubChem lookup (network, on request)");
    a.dict->require_subcommand(0, 1);
    a.dictList = a.dict->add_subcommand("list", "residue codes with their SMILES");
    a.dictList->add_flag("--json", a.listJson, "machine-readable output");
    a.dictSmiles = a.dict->add_subcommand("smiles", "print the SMILES given to acedrg / eLBOW for a code");
    a.dictSmiles->add_option("code", a.code, "residue code")->required();
    a.dictFinish = a.dict->add_subcommand("finish", "rename an acedrg / eLBOW dictionary to the model's atom names");
    a.dictFinish->add_option("input", a.finishInputs, "dictionary .cif, or folders / glob patterns to search")
        ->required();
    lastWins(a.dictFinish->add_option("--code", a.code, "residue code, e.g. 0A5"))->required();
    lastWins(a.dictFinish->add_option("--output,-o", a.finishOutput, "renamed dictionary"))->required();
    lastWins(a.dictFinish->add_option("--group", a.group,
                                      "_chem_comp.group (default RNA: polymer linking in Coot / Phenix; 'non-polymer' "
                                      "for a free ligand)"));
    lastWins(a.dictFinish->add_option("--name", a.name, "_chem_comp.name (default: from the table)"));
    lastWins(a.dictFinish->add_option("--source-id", a.sourceId, "component id in the input when it holds several"));
    a.dictFinish->add_flag("--allow-protonation", a.allowProtonation,
                           "accept hydrogen counts that differ from the reference outside the phosphate");
    a.dictFinish->add_flag("--keep-names", a.keepNames,
                           "only set id / name / group, keep the atom names (codes without reference names)");
    lastWins(a.dictFinish->add_option("--project,-p,--params", a.finishProject,
                                      "project file whose step log records the step"));
    a.dictPubchem = a.dict->add_subcommand(
        "pubchem", "fetch a SMILES from PubChem (network) and check it against the known codes");
    a.dictPubchem->add_option("compound", a.compound, "PubChem CID or compound name")->required();
    lastWins(a.dictPubchem->add_option("--code", a.pubchemCode, "only check against this code"));
    floatOption(*a.dictPubchem, "--timeout", a.timeout, "seconds (default 15)");

    // ---- density fit / simulate (cli.py)
    a.densityFit = density.add_subcommand("fit", "fit the structure into the project's (or --map) density map");
    addCommon(*a.densityFit, a.fitCommon);
    lastWins(a.densityFit->add_option("--map", a.fitMap, "MRC / CCP4 map (default: density.map_path of the project)"));
    floatOption(*a.densityFit, "--resolution", a.resolution, "map resolution in Å (simulated-map width)");
    floatOption(*a.densityFit, "--level", a.level, "contour level in σ (shadow and centring)");
    a.frameOption = lastWins(
        a.densityFit->add_option("--frame", a.frame, "'ref' (first analysed frame) or a trajectory frame number"));
    choice(lastWins(a.densityFit->add_option(
               "--search", a.search, "global search, or refine the stored placement (default: density.search)")),
           {"global", "manual"});
    floatOption(*a.densityFit, "--step", a.step, "angular step of the global search (degrees)");
    a.saveOption = lastWins(a.densityFit->add_option("--save", a.save,
                                                     "store the placement in the project (--project, or the given path)"))
                       ->expected(0, 1)
                       ->type_name("PRESET");
    a.densitySimulate = density.add_subcommand("simulate", "write a map simulated from one frame (molmap-style Gaussians)");
    addCommon(*a.densitySimulate, a.simulateCommon);
    floatOption(*a.densitySimulate, "--resolution", a.resolution, "resolution in Å")->required();
    floatOption(*a.densitySimulate, "--spacing", a.spacing, "grid spacing in Å (default resolution / 3)");
    lastWins(a.densitySimulate->add_option("--frame", a.simulateFrame, "'ref' or a trajectory frame number"));
    a.densitySimulate->add_flag("--placed", a.placed, "apply the stored placement (map frame)");
    lastWins(a.densitySimulate->add_option("--output,-O", a.simulateOutput, "output .mrc"))->required();
}

ChemCommands::~ChemCommands() = default;

bool ChemCommands::parsed() const {
    const Arguments& a = *m_args;
    return a.pdb->parsed() || a.build->parsed() || a.dict->parsed() || densityParsed();
}

bool ChemCommands::densityParsed() const { return m_args->densityFit->parsed() || m_args->densitySimulate->parsed(); }

int ChemCommands::run(Console& console, const std::vector<std::string>& arguments) {
    Arguments& a = *m_args;
    if (a.densityFit->parsed()) {
        return densityFitRun(a, *a.densityFit, console);
    }
    if (a.densitySimulate->parsed()) {
        return densitySimulateRun(a, *a.densitySimulate, console);
    }
    for (CLI::App* group : {a.pdb, a.build, a.dict}) {
        if (group->parsed() && group->get_subcommands().empty()) {
            console.out << group->help();
            return 2;
        }
    }
    if (a.pdb->parsed()) {
        return pdbRun(a, console, arguments);
    }
    if (a.build->parsed()) {
        return buildRun(a, console, arguments);
    }
    return dictRun(a, console, arguments);
}

}  // namespace mdbridge::app::detail
