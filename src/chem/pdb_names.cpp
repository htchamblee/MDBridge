// Internal adapter (pdb_names.hpp): MDATK's pdbtools.convert_names / set_records / scheme_name for
// the helix builder, ported from mdatk/pdbtools.py (v0.0.4.1; _nucleotide with CH-05 of 0.1.0).

#include "pdb_names.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/json.hpp"

namespace mdbridge::chem::detail {

namespace {

using io::PdbResidue;
using io::PdbStructure;

constexpr std::array<std::string_view, 4> kSchemeNames = {"amber", "ccd", "rosetta", "pdb"};
constexpr std::array<std::string_view, 3> kRecordPolicies = {"keep", "rosetta", "pdb"};

constexpr std::array<std::string_view, 20> kAmino = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                                     "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                                     "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
constexpr std::array<std::string_view, 8> kAminoVariants = {"HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN"};
constexpr std::array<std::string_view, 13> kStandardExtra = {"A",  "C",  "G",  "U",  "I",  "N",  "DA",
                                                             "DC", "DG", "DT", "DI", "DN", "UNK"};

template <std::size_t N>
bool contains(const std::array<std::string_view, N>& items, std::string_view value) {
    return std::find(items.begin(), items.end(), value) != items.end();
}

std::string asciiUpper(std::string_view text) {
    std::string out(text);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return out;
}

std::string primed(std::string_view name) {
    std::string out(name);
    std::replace(out.begin(), out.end(), '*', '\'');
    return out;
}

// pdbtools._nucleotide, with CH-05 (MDATK 0.1.0): an RNA name on a 2'-deoxy sugar is DNA
std::optional<io::ResidueName> nucleotide(const PdbResidue& res) {
    std::optional<io::ResidueName> info = io::parseResidueName(res.name);
    if (!info) {
        return std::nullopt;
    }
    bool backbone = false;
    for (const io::PdbAtom& atom : res.atoms) {
        if (atom.name.find('\'') != std::string::npos || atom.name.find('*') != std::string::npos ||
            atom.name == "P") {
            backbone = true;
            break;
        }
    }
    if (!backbone) {
        return std::nullopt;
    }
    if (info->sugar == "rna") {
        bool c2 = false;
        bool o2 = false;
        for (const io::PdbAtom& atom : res.atoms) {
            const std::string name = primed(atom.name);
            c2 = c2 || name == "C2'";
            o2 = o2 || name == "O2'";
        }
        if (c2 && !o2) {
            info->sugar = "dna";
        }
    }
    return info;
}

bool isPolymer(const PdbResidue& res) {
    const std::string upper = asciiUpper(res.name);
    return nucleotide(res).has_value() || contains(kAmino, upper) || contains(kAminoVariants, upper);
}

bool isStandardName(std::string_view name) {
    return contains(kAmino, name) || contains(kStandardExtra, name);
}

bool hasAtom(const PdbResidue& res, std::string_view name) {
    return std::any_of(res.atoms.begin(), res.atoms.end(), [&](const io::PdbAtom& a) { return a.name == name; });
}

// pdbtools._listing
std::string listing(const std::vector<std::string>& items, std::size_t limit = 8) {
    std::string out;
    for (std::size_t k = 0; k < items.size() && k < limit; ++k) {
        if (k) {
            out += ", ";
        }
        out += items[k];
    }
    if (items.size() > limit) {
        out += " and " + std::to_string(items.size() - limit) + " more";
    }
    return out;
}

template <std::size_t N>
std::string unknown(std::string_view what, std::string_view value, const std::array<std::string_view, N>& names) {
    std::string out = "unknown " + std::string(what) + " " + core::pythonStringRepr(std::string(value)) + " (choose from ";
    for (std::size_t k = 0; k < N; ++k) {
        out += (k ? ", " : "") + std::string(names[k]);
    }
    return out + ")";
}

// A Counter: insertion order kept
struct Counts {
    std::vector<std::pair<std::pair<std::string, std::string>, int>> items;
    void add(const std::string& a, const std::string& b) {
        for (auto& item : items) {
            if (item.first.first == a && item.first.second == b) {
                ++item.second;
                return;
            }
        }
        items.push_back({{a, b}, 1});
    }
};

}  // namespace

std::string schemeName(char base, char chirality, std::string_view terminal, std::string_view scheme) {
    const std::string b(1, base);
    if (scheme == "amber") {
        if (chirality == 'L') {  // nucleic12_LRNA.lib defines LX, LX5 and LX3 only
            return "L" + b + (terminal == "5" || terminal == "N" ? "5" : terminal == "3" ? "3" : "");
        }
        return b + std::string(terminal);
    }
    if (scheme == "ccd") {
        if (chirality == 'L') {
            return "0" + b + (terminal == "5" ? "5" : terminal == "3" ? "3" : "");
        }
        return b;
    }
    if (scheme == "rosetta" || scheme == "pdb") {
        return chirality == 'L' ? "0" + b : b;
    }
    throw io::PdbError(unknown("residue name scheme", scheme, kSchemeNames));
}

std::vector<std::string> convertNames(PdbStructure& structure, std::string_view scheme, bool terminals) {
    if (!contains(kSchemeNames, scheme)) {
        throw io::PdbError(unknown("residue name scheme", scheme, kSchemeNames));
    }
    Counts changes;
    std::vector<std::string> termMismatch, phosphate, singleL, keptPhosphate;
    for (std::vector<PdbResidue>& seg : structure.segments) {
        // _terminals: 5' / 3' / N over the nucleotides of the segment
        std::vector<std::size_t> nts;
        for (std::size_t k = 0; k < seg.size(); ++k) {
            if (nucleotide(seg[k])) {
                nts.push_back(k);
            }
        }
        std::vector<std::string> terms(seg.size());
        if (nts.size() == 1) {
            terms[nts[0]] = "N";
        } else if (!nts.empty()) {
            terms[nts.front()] = "5";
            terms[nts.back()] = "3";
        }
        for (std::size_t k = 0; k < seg.size(); ++k) {
            PdbResidue& res = seg[k];
            const std::optional<io::ResidueName> nt = nucleotide(res);
            if (!nt || nt->sugar != "rna") {
                continue;
            }
            const std::string segTerm = terminals ? terms[k] : "";
            std::string term = segTerm;
            const char chir = res.chirality && !res.chirality->empty() ? (*res.chirality)[0] : nt->chirality;
            if (scheme == "ccd" && chir == 'L' && (term == "5" || term == "N") && hasAtom(res, "P")) {
                term.clear();  // a 5' end with its phosphate is the CCD component itself
                keptPhosphate.push_back(res.label());
            }
            const std::string name = schemeName(nt->base, chir, term, scheme);
            const bool singleLUnit = segTerm == "N" && chir == 'L' && nt->terminal == "5";
            if (terminals && !nt->terminal.empty() && nt->terminal != segTerm && !singleLUnit) {
                termMismatch.push_back(res.label());
            }
            if (name != res.name) {
                changes.add(res.name, name);
                res.name = name;
            }
            if (scheme == "amber" && (terms[k] == "5" || terms[k] == "N") && hasAtom(res, "P")) {
                phosphate.push_back(res.label());
            }
            if (scheme == "amber" && term == "N" && chir == 'L') {
                singleL.push_back(res.label());
            }
        }
    }
    std::vector<std::string> notes;
    if (!changes.items.empty()) {
        auto items = changes.items;
        std::stable_sort(items.begin(), items.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        int total = 0;
        std::string summary;
        for (const auto& [key, n] : items) {
            total += n;
            summary += (summary.empty() ? "" : ", ") + key.first + "->" + key.second + " x" + std::to_string(n);
        }
        notes.push_back("Renamed " + std::to_string(total) + " residues to the " + std::string(scheme) +
                        " scheme: " + summary);
    } else {
        notes.push_back("Residue names already follow the " + std::string(scheme) + " scheme");
    }
    if (!termMismatch.empty()) {
        notes.push_back("Terminal names that are not at a segment end (missing TER?): " + listing(termMismatch));
    }
    if (!phosphate.empty()) {
        notes.push_back("WARNING: 5'-terminal residues carry P/OP1/OP2, which the Amber 5' units do not have "
                        "(tleap cannot type them; remove them with --strip-5p or cap the end): " +
                        listing(phosphate));
    }
    if (!keptPhosphate.empty()) {
        notes.push_back("5'-terminal L residues with a 5'-phosphate keep the CCD component names (0A ..., which "
                        "include P/OP1/OP2/OP3; the 5' variants 0A5 ... have a 5'-OH): " +
                        listing(keptPhosphate));
    }
    if (!singleL.empty()) {
        notes.push_back("nucleic12_LRNA.lib has no unit for a single-residue L segment; wrote the 5' unit "
                        "(its O3' lacks a hydrogen): " +
                        listing(singleL));
    }
    return notes;
}

std::vector<std::string> setRecords(PdbStructure& structure, std::string_view policy) {
    if (!contains(kRecordPolicies, policy)) {
        throw io::PdbError(unknown("record policy", policy, kRecordPolicies));
    }
    if (policy == "keep") {
        return {};
    }
    Counts changed;
    for (std::vector<PdbResidue>& seg : structure.segments) {
        for (PdbResidue& res : seg) {
            std::string record;
            if (policy == "rosetta") {
                record = isPolymer(res) ? "ATOM" : "";
            } else {
                record = isStandardName(res.name) ? "ATOM" : "HETATM";
            }
            if (record.empty() || std::all_of(res.atoms.begin(), res.atoms.end(),
                                              [&](const io::PdbAtom& a) { return a.record == record; })) {
                continue;
            }
            changed.add(res.name, record);
            for (io::PdbAtom& atom : res.atoms) {
                atom.record = record;
            }
        }
    }
    if (changed.items.empty()) {
        return {};
    }
    auto items = changed.items;
    std::sort(items.begin(), items.end());
    std::string summary;
    for (const auto& [key, n] : items) {
        summary += (summary.empty() ? "" : ", ") + key.first + "->" + key.second + " x" + std::to_string(n);
    }
    return {"Record policy " + std::string(policy) + ": " + summary};
}

}  // namespace mdbridge::chem::detail
