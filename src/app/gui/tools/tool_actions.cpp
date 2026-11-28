#include "gui/tools/tool_actions.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <map>
#include <system_error>
#include <tuple>
#include <utility>

#include "commands.hpp"
#include "mdbridge/chem/fit_to_map.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/mutate.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/step_log.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::app::tools {

namespace fs = std::filesystem;
using core::Value;

namespace {

std::string joined(const std::vector<std::string>& lines, const char* separator = "\n") {
    std::string out;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        out += (i ? separator : "") + lines[i];
    }
    return out;
}

// open(path, "w", encoding="utf-8").write(text), as the command line writes (chem_commands.cpp):
// OSError texts name the path as given
void writeText(const std::string& path, const std::string& text) {
    const fs::path p = core::pathFromUtf8(path);
    std::ofstream out(p, std::ios::binary | std::ios::trunc);
    if (!out) {
        const int error = errno != 0 ? errno : EIO;
        std::error_code code(error, std::generic_category());
        std::error_code ec;
        if (fs::is_directory(p, ec)) {
            code = std::make_error_code(std::errc::is_a_directory);
        }
        throw fs::filesystem_error("open", p, code);
    }
    out << text;
    out.close();
    if (!out) {
        throw fs::filesystem_error("write", p, std::make_error_code(std::errc::io_error));
    }
}

std::string absolute(const std::string& path) {
#ifdef _WIN32
    return core::pathToUtf8(fs::absolute(core::pathFromUtf8(path)).lexically_normal());
#else
    return core::posixpath::abspath(path, core::currentDirectory());
#endif
}

bool isSeparator(char c) {
#ifdef _WIN32
    return c == '/' || c == '\\';
#else
    return c == '/';
#endif
}

std::vector<std::string> keysFor(std::initializer_list<std::string> resources) {
    const std::vector<std::string> ids(resources);
    return core::citations::keysFor(ids);
}

const Value* member(const Value& table, std::string_view key) {
    const Value* v = table.find(key);
    return v;
}

std::int64_t integerOf(const Value* v) {
    const auto* i = v ? v->get<std::int64_t>() : nullptr;
    return i ? *i : 0;
}

std::string stringOf(const Value* v) { return v ? v->stringOr() : std::string(); }

}  // namespace

std::string errorText(const std::exception& error) { return pythonErrorText(error); }

bool isFile(const std::string& path) {
    if (path.empty()) {
        return false;
    }
    std::error_code ec;
    return fs::is_regular_file(core::pathFromUtf8(path), ec);
}

bool samePath(const std::string& a, const std::string& b) { return absolute(a) == absolute(b); }

std::string suggestedOutput(const std::string& input, std::string_view suffix) {
    if (input.empty()) {
        return {};
    }
    // genericpath._splitext: the last dot of the file name, unless the name is only leading dots
    std::size_t nameStart = 0;
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (isSeparator(input[i])) nameStart = i + 1;
    }
#ifdef _WIN32
    if (nameStart == 0 && input.size() >= 2 && input[1] == ':') nameStart = 2;  // ntpath's drive
#endif
    std::string stem = input;
    std::string ext;
    const std::size_t dot = input.rfind('.');
    if (dot != std::string::npos && dot >= nameStart) {
        for (std::size_t k = nameStart; k < dot; ++k) {
            if (input[k] != '.') {
                stem = input.substr(0, dot);
                ext = input.substr(dot);
                break;
            }
        }
    }
    return stem + "_" + std::string(suffix) + (ext.empty() ? std::string(".pdb") : ext);
}

std::optional<std::vector<std::string>> splitIds(const std::string& text) {
    std::vector<std::string> out;
    std::size_t begin = 0;
    while (true) {
        const std::size_t comma = text.find(',', begin);
        const std::string part = core::pythonStrip(
            std::string_view(text).substr(begin, comma == std::string::npos ? std::string::npos : comma - begin));
        if (!part.empty()) out.push_back(part);
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    if (out.empty()) {
        return std::nullopt;
    }
    return out;
}

// ---- PDB clean-up ------------------------------------------------------------------------------

std::string pdbCheckReport(const std::string& input, const std::string& lRanges) {
    const Value info = chem::check(chem::PdbSource::file(core::pathFromUtf8(input)),
                                   lRanges.empty() ? chem::RangeSpec() : chem::RangeSpec(lRanges));
    std::vector<std::string> lines;
    lines.push_back(core::pathBasename(input) + ": " + std::to_string(integerOf(member(info, "segments"))) +
                    " TER segment(s), " + std::to_string(integerOf(member(info, "conect_records"))) +
                    " CONECT record(s)");
    lines.emplace_back();
    if (const auto* chains = info.find("chains") ? info.find("chains")->get<core::Table>() : nullptr) {
        for (const auto& [chain, c] : *chains) {
            std::vector<std::pair<std::string, std::int64_t>> calls;
            if (const auto* chirality = c.find("chirality") ? c.find("chirality")->get<core::Table>() : nullptr) {
                for (const auto& [call, n] : *chirality) calls.emplace_back(call, integerOf(&n));
            }
            std::sort(calls.begin(), calls.end());
            std::string chir;
            for (std::size_t k = 0; k < calls.size(); ++k) {
                chir += (k ? ", " : "") + std::to_string(calls[k].second) + " " + calls[k].first;
            }
            std::vector<std::string> names;
            if (const auto* list = c.find("names") ? c.find("names")->get<core::Array>() : nullptr) {
                for (const Value& n : *list) names.push_back(n.stringOr());
            }
            lines.push_back("Chain " + chain + ": " + std::to_string(integerOf(c.find("residues"))) + " residues (" +
                            stringOf(c.find("numbers")) + "), names " + joined(names, " ") + "; chirality " +
                            (chir.empty() ? std::string("-") : chir));
            lines.push_back("  " + stringOf(c.find("sequence")));
        }
    }
    if (const Value* other = info.find("other_residues")) {
        if (const auto* table = other->get<core::Table>(); table && !table->empty()) {
            lines.push_back("Other residues: " + core::toJsonAscii(*other, -1));  // json.dumps
        }
    }
    lines.emplace_back();
    if (const auto* notes = info.find("notes") ? info.find("notes")->get<core::Array>() : nullptr) {
        for (const Value& n : *notes) lines.push_back(n.stringOr());
    }
    return joined(lines);
}

ToolResult cleanPdb(const PdbCleanRequest& r, const fs::path& project) {
    ToolResult result;
    result.output = r.output;
    core::logged(project, "pdb clean", "PDB clean-up of " + core::pathBasename(r.input), {}, [&](core::LogEntry& entry) {
        chem::CleanOptions o;
        if (!r.scheme.empty()) o.scheme = r.scheme;
        if (!r.records.empty()) o.records = r.records;
        o.stripConect = r.stripConect;
        o.detectL = r.detectL;
        o.lRanges = r.lRanges.empty() ? chem::RangeSpec() : chem::RangeSpec(r.lRanges);
        if (!r.chainReference.empty()) o.chainReference = chem::PdbSource::file(core::pathFromUtf8(r.chainReference));
        o.chainIds = r.chainIds;
        o.assignChains = r.assignChains;
        o.copyNumbers = r.copyNumbers;
        o.splitSegments = r.splitSegments;
        o.splitAt = r.splitAt.empty() ? chem::RangeSpec() : chem::RangeSpec(r.splitAt);
        o.combine = r.combine;
        o.combineInto = r.combineInto;
        o.sortChains = r.sortChains;
        o.renumber = r.renumber;
        o.start = r.start;
        o.perChain = r.perChain;
        o.keepGaps = r.keepGaps;
        o.seqres = r.seqres;
        o.makeFasta = !r.fasta.empty();
        o.terminalNames = r.terminalNames;
        o.strip5p = r.strip5p;
        const chem::CleanResult res = chem::clean(chem::PdbSource::file(core::pathFromUtf8(r.input)), o);
        std::string text = res.text;
        std::vector<std::string> notes = res.notes;
        if (r.annotation) {
            const chem::AnnotationResult ann = chem::annotate(chem::PdbSource::text(text), *r.annotation);
            text = ann.text;
            for (const std::string& n : ann.notes) {
                if (std::find(notes.begin(), notes.end(), n) == notes.end()) notes.push_back(n);
            }
            // CH-14: never Barnaba, so never its citation
        }
        writeText(r.output, text);
        entry.files = {r.output};
        if (!r.fasta.empty() && res.fasta && !res.fasta->empty()) {
            writeText(r.fasta, *res.fasta);
            entry.files.push_back(r.fasta);
            result.fasta = r.fasta;
        }
        entry.summary = joined(notes);
        result.notes = std::move(notes);
    });
    return result;
}

// ---- builder -----------------------------------------------------------------------------------

namespace {

std::string joinedChains(const std::vector<std::string>& chains) {
    // MDATK passes the chain ids as a tuple: an id typed with a comma is one id, longer than one
    // character ("chain ids must be single characters"), not two
    std::vector<std::string> parts;
    for (std::string c : chains) {
        for (std::size_t at = c.find(','); at != std::string::npos; at = c.find(',', at + 2)) c.replace(at, 1, "??");
        parts.push_back(std::move(c));
    }
    return joined(parts, ",");
}

std::string joinedNumbers(const std::vector<std::int64_t>& numbers) {
    std::vector<std::string> parts;
    for (const std::int64_t n : numbers) parts.push_back(std::to_string(n));
    return joined(parts, ",");
}

std::optional<std::string_view> optionalText(const std::string& text) {
    return text.empty() ? std::nullopt : std::optional<std::string_view>(text);
}

}  // namespace

bool helixInBackground(const HelixRequest& request) {
    return request.pairing != "wc" || !(request.triplex.empty() || request.triplex == "none");
}

std::string helixCheckReport(const HelixRequest& r) {
    std::vector<std::string> lines;
    for (const chem::HelixMessage& m : chem::checkHelix(r.sequence, r.partner, r.pairing, optionalText(r.triplex), r.third)) {
        std::string level = m.level;
        std::transform(level.begin(), level.end(), level.begin(), [](char c) {
            return c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : c;
        });
        lines.push_back(level + ": " + m.text);
    }
    return joined(lines);
}

ToolResult buildHelix(const HelixRequest& r, const fs::path& project, const core::Reporter& reporter) {
    reporter.checkCancel();
    ToolResult result;
    result.output = r.output;
    const std::string title = "Ideal " + chem::helixKind(r.pairing, optionalText(r.triplex)) + " " + r.chirality +
                              "-RNA helix " + core::pythonStrip(r.sequence);
    core::logged(project, "build helix", title, {}, [&](core::LogEntry& entry) {
        chem::HelixOptions o;
        o.partner = r.partner;
        o.chirality = r.chirality;
        o.chains = joinedChains(r.chains);
        o.start = joinedNumbers(r.start);
        if (!r.scheme.empty()) o.scheme = r.scheme;
        o.terminals = r.terminals;
        o.hydrogens = r.hydrogens;
        o.phosphate5 = r.phosphate5;
        o.pairing = r.pairing;
        if (!r.triplex.empty()) o.triplex = r.triplex;
        o.third = r.third;
        std::vector<std::string> notes;
        o.reporter = &reporter;  // the fit and the closure stop on the job's cancel
        const std::string text = chem::buildHelix(r.sequence, o, &notes);
        reporter.checkCancel();  // a build abandoned at its end writes nothing and logs no step
        writeText(r.output, text);
        entry.files = {r.output};
        entry.summary = joined(notes);
        // D19: the Hoogsteen duplex and the pyrimidine-motif triplex also cite their RNA triple helices
        const bool triplexes = r.pairing == "hoogsteen" || (r.pairing == "wc" && r.triplex == "hoogsteen");
        entry.citations = triplexes                ? keysFor({"arnott-arna", "pdb-4jrd", "pdb-triplexes"})
                          : helixInBackground(r) ? keysFor({"arnott-arna", "pdb-4jrd"})
                                                 : keysFor({"arnott-arna"});
        result.notes = std::move(notes);
    });
    return result;
}

std::string mutateCheckReport(const std::string& input, const std::string& mutations) {
    const std::vector<chem::Mutation> muts = chem::parseMutations(mutations);
    const io::PdbStructure s = chem::PdbSource::file(core::pathFromUtf8(input)).read();
    std::map<std::tuple<std::string, std::int64_t, std::string>, const io::PdbResidue*> have;
    for (const io::PdbResidue* res : s.residues()) {
        have[{res->chain, res->resseq, res->icode}] = res;  // the last of equal ids, as MDATK's dict
    }
    std::vector<std::string> lines;
    for (const chem::Mutation& m : muts) {
        const auto it = have.find({m.chain, m.resseq, m.icode});
        lines.push_back(m.where() + " -> " + std::string(1, m.base) + ": " +
                        (it != have.end() ? it->second->name : std::string("NOT FOUND")));
    }
    return joined(lines);
}

ToolResult mutate(const std::string& input, const std::string& mutations, bool compensatory, bool wobble,
                  const std::string& output, const fs::path& project) {
    ToolResult result;
    result.output = output;
    const std::string spec = core::pythonStrip(mutations);
    core::logged(project, "build mutate", "Mutation of " + core::pathBasename(input) + ": " + spec, {},
                 [&](core::LogEntry& entry) {
                     chem::MutateOptions o;
                     o.compensatory = compensatory;
                     o.wobble = wobble;
                     chem::MutateResult res = chem::mutate(chem::PdbSource::file(core::pathFromUtf8(input)), spec, o);
                     writeText(output, res.text);
                     entry.files = {output};
                     entry.summary = joined(res.notes);
                     result.notes = std::move(res.notes);
                 });
    return result;
}

ToolResult fitToMap(const std::string& input, const std::string& map, double resolution, const std::string& search,
                    double step, const std::string& output, const fs::path& project, const core::Reporter& reporter) {
    ToolResult result;
    result.output = output;
    const std::string title = "Rigid-body fit of " + core::pathBasename(input) + " into " + core::pathBasename(map);
    core::logged(project, "build fit", title, {}, [&](core::LogEntry& entry) {
        chem::FitToMapOptions o;
        o.search = search;
        o.angularStep = step;
        o.control.reporter = &reporter;
        chem::FitToMapResult res =
            chem::fitToMap(chem::PdbSource::file(core::pathFromUtf8(input)), core::pathFromUtf8(map), resolution, o);
        reporter.checkCancel();  // Cancel pressed as the fit ended: nothing is written (MDATK 0.0.13)
        writeText(output, res.text);
        entry.files = {output};
        entry.summary = joined(res.notes);
        result.notes = std::move(res.notes);
        result.cc = res.cc;
    });
    return result;
}

}  // namespace mdbridge::app::tools
