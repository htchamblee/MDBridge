// MDATK v0.0.4.1's base-pair annotation (mdatk/annotate.py, MDATK's own code) on io's PDB model and
// the builtin detector of secstruct.cpp, with D13 and CH-14 (annotate.hpp). Every loop follows
// MDATK's order, so records, notes and text are MDATK's byte for byte (tests/chem/test_annotate.cpp).

#include "mdbridge/chem/annotate.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::chem {

using io::PdbResidue;
using io::PdbStructure;

namespace {

constexpr std::string_view kHeaderRecord = "REMARK 950 MDBridge BASE-PAIR ANNOTATION";  // D13
constexpr std::string_view kStackingNote = "The builtin engine finds canonical base pairs only (no stacking)";

bool startsWith(std::string_view text, std::string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

// f"{x:g}" (Python's general format, precision 6)
std::string formatGeneral(double x) {
    if (std::isnan(x)) {
        return "nan";
    }
    char text[32];
    std::snprintf(text, sizeof text, "%g", x);
    return text;
}

// (s or fill)[:1]
std::u32string firstCharacter(const std::string& text, char32_t fill) {
    const std::u32string chars = core::decodeUtf8(text);
    return chars.empty() ? std::u32string(1, fill) : chars.substr(0, 1);
}

void checkEngine(const AnnotateOptions& options) {
    if (std::find(kAnnotationEngines.begin(), kAnnotationEngines.end(), options.engine) ==
        kAnnotationEngines.end()) {
        throw core::ValueError("unknown annotation engine " + core::pythonStringRepr(options.engine) +
                               " (choose from auto, builtin, barnaba)");
    }
}

std::string pairCode(int kind) {
    return kind == kWatsonCrick ? "WCc" : "GUc";
}

std::vector<std::string> dotbracketRecords(const std::vector<const PdbResidue*>& residues,
                                           const std::vector<std::vector<std::size_t>>& strands,
                                           const std::vector<ResiduePair>& pairs) {
    std::vector<std::string> out;
    for (const std::vector<std::size_t>& group : annotationStrandGroups(strands, pairs)) {
        std::vector<std::size_t> rows;
        std::vector<std::int64_t> breaks;
        for (const std::size_t c : group) {
            if (!rows.empty()) {
                breaks.push_back(static_cast<std::int64_t>(rows.size()));
            }
            rows.insert(rows.end(), strands[c].begin(), strands[c].end());
        }
        std::map<std::size_t, std::int64_t> local;
        for (std::size_t n = 0; n < rows.size(); ++n) {
            local.emplace(rows[n], static_cast<std::int64_t>(n));
        }
        std::vector<ResiduePair> groupPairs;
        for (const auto& [i, j] : pairs) {
            const auto a = local.find(static_cast<std::size_t>(i));
            const auto b = local.find(static_cast<std::size_t>(j));
            if (a != local.end() && b != local.end()) {
                groupPairs.emplace_back(a->second, b->second);
            }
        }
        std::string seq;
        for (std::size_t n = 0; n < rows.size(); ++n) {
            if (std::find(breaks.begin(), breaks.end(), static_cast<std::int64_t>(n)) != breaks.end()) {
                seq += kChainSeparator;
            }
            seq += baseLetter(residues[rows[n]]->name).value_or('?');  // a nucleotide always has one
        }
        std::vector<std::string> chains;  // dict.fromkeys: unique, in order
        for (const std::size_t c : group) {
            const std::string& id = residues[strands[c][0]]->chain;
            std::string chain = id.empty() ? std::string("_") : id;
            if (std::find(chains.begin(), chains.end(), chain) == chains.end()) {
                chains.push_back(std::move(chain));
            }
        }
        std::string joined;
        for (std::size_t k = 0; k < chains.size(); ++k) {
            joined += (k ? "," : "") + chains[k];
        }
        out.push_back(std::string(kAnnotationRemark) + " CHAINS " + joined);
        for (std::string& line : annotationWrapped("SEQ", seq)) {
            out.push_back(std::move(line));
        }
        const std::string dbn = pairsToDotbracket(groupPairs, static_cast<std::int64_t>(rows.size()), breaks);
        for (std::string& line : annotationWrapped("DBN", dbn)) {
            out.push_back(std::move(line));
        }
    }
    return out;
}

Notes gapNotes(const std::vector<const PdbResidue*>& residues, const std::vector<std::vector<std::size_t>>& strands) {
    std::vector<std::string> gaps;
    for (const std::vector<std::size_t>& rows : strands) {
        for (std::size_t k = 0; k + 1 < rows.size(); ++k) {
            const PdbResidue& ra = *residues[rows[k]];
            const PdbResidue& rb = *residues[rows[k + 1]];
            if (rb.resseq - ra.resseq > 1) {
                gaps.push_back((ra.chain.empty() ? std::string("_") : ra.chain) + ":" + std::to_string(ra.resseq) +
                               "-" + std::to_string(rb.resseq));
            }
        }
    }
    if (gaps.empty()) {
        return {};
    }
    std::string joined;
    for (std::size_t k = 0; k < gaps.size(); ++k) {
        joined += (k ? ", " : "") + gaps[k];
    }
    return {"Residue numbering gaps (" + joined + "): missing residues are not in SEQ / DBN, so the residues on "
            "either side are written next to each other"};
}

// annotate.py:145-231 after the reading: `notes` holds the read notes as annotate reports them
struct Annotated {
    Notes notes;
    std::vector<std::string> records;
};

Annotated annotateStructure(const PdbStructure& s, Notes notes, const AnnotateOptions& o) {
    const std::vector<const PdbResidue*> residues = s.residues();
    const std::vector<std::vector<std::size_t>> strands = annotationStrands(s);
    std::set<std::size_t> nucleotides;
    for (const std::vector<std::size_t>& rows : strands) {
        nucleotides.insert(rows.begin(), rows.end());
    }
    std::vector<double> xyz;
    std::vector<AtomIndex> atomIndex;
    std::vector<std::string> names;  // '' never pairs
    atomIndex.reserve(residues.size());
    names.reserve(residues.size());
    std::int64_t nAtoms = 0;
    for (std::size_t k = 0; k < residues.size(); ++k) {
        AtomIndex index;
        for (const auto& [name, p] : residues[k]->coords()) {
            index.emplace(name, nAtoms++);
            xyz.insert(xyz.end(), p.begin(), p.end());
        }
        atomIndex.push_back(std::move(index));
        names.push_back(nucleotides.count(k) ? residues[k]->name : std::string());
    }
    if (o.engine == "barnaba") {  // CH-14: before the engine's notes
        notes.emplace_back(kBarnabaNote);
    }
    const PairDetector detector(names, atomIndex, o.cutoff, o.wobble);
    const std::vector<BasePair> found = detector.detect(xyz);
    if (o.stacking) {
        notes.emplace_back(kStackingNote);
    }
    std::size_t watsonCrick = 0;
    std::size_t wobble = 0;
    std::size_t other = 0;
    for (const BasePair& p : found) {
        const int kind = lwKind(pairCode(p.kind));
        (kind == kWatsonCrick ? watsonCrick : kind == kWobble ? wobble : other) += 1;
    }
    const std::string remark(kAnnotationRemark);
    std::vector<std::string> records = {std::string(kHeaderRecord),
                                        remark + " ENGINE builtin, canonical pairs from H-bond distances <= " +
                                            formatGeneral(o.cutoff) + " A"};
    records.push_back(remark + " PAIRS " + std::to_string(found.size()) + " (WATSON-CRICK " +
                      std::to_string(watsonCrick) + ", WOBBLE " + std::to_string(wobble) + ", OTHER " +
                      std::to_string(other) + ")");
    if (o.pairs) {
        for (const BasePair& p : found) {
            records.push_back(remark + " PAIR  " + annotationResidueField(*residues[static_cast<std::size_t>(p.i)]) +
                              " " + annotationResidueField(*residues[static_cast<std::size_t>(p.j)]) + " " +
                              pairCode(p.kind));
        }
    }
    if (o.dotbracket) {
        std::vector<ResiduePair> canonical;
        for (const BasePair& p : found) {
            const int kind = lwKind(pairCode(p.kind));
            if (kind == kWatsonCrick || kind == kWobble) {
                canonical.emplace_back(p.i, p.j);
            }
        }
        UniquePairs unique = uniquePairs(canonical);
        if (!unique.dropped.empty()) {
            notes.push_back("Left " + std::to_string(unique.dropped.size()) +
                            " canonical pair(s) out of the dot-bracket: their residues already pair with another base");
        }
        for (std::string& line : dotbracketRecords(residues, strands, unique.keep)) {
            records.push_back(std::move(line));
        }
        for (std::string& note : gapNotes(residues, strands)) {
            notes.push_back(std::move(note));
        }
    }
    notes.push_back("Annotated " + std::to_string(found.size()) + " base pair(s) with the builtin engine");
    return {std::move(notes), std::move(records)};
}

// "Kept model 1 of N" as annotate reports it
std::string modelNote(const std::string& note) {
    return "Annotated model 1 of " + note.substr(note.rfind(' ') + 1) + " (all models are copied unchanged)";
}

}  // namespace

std::string annotationResidueField(const PdbResidue& residue) {
    std::u32string out = firstCharacter(residue.chain, U' ');
    std::string number = std::to_string(residue.resseq);
    if (number.size() < 4) {
        number.insert(0, 4 - number.size(), ' ');
    }
    out += core::decodeUtf8(number);
    out += firstCharacter(residue.icode, U' ');
    std::u32string name = core::decodeUtf8(residue.name);
    if (name.size() < 3) {
        name.append(3 - name.size(), U' ');
    }
    out += name;
    return core::encodeUtf8(out);
}

std::vector<std::string> annotationWrapped(std::string_view keyword, std::string_view text) {
    const std::u32string chars = core::decodeUtf8(text);
    std::string key(keyword);
    if (core::characterCount(key) < 5) {
        key.append(5 - core::characterCount(key), ' ');
    }
    std::vector<std::string> out;
    for (std::size_t k = 0; k < std::max<std::size_t>(1, chars.size()); k += kAnnotationWrap) {
        out.push_back(std::string(kAnnotationRemark) + " " + key + " " +
                      core::encodeUtf8(chars.substr(std::min(k, chars.size()), kAnnotationWrap)));
    }
    return out;
}

std::vector<std::vector<std::size_t>> annotationStrands(const PdbStructure& structure) {
    std::vector<std::vector<std::size_t>> out;
    std::size_t k = 0;
    for (const std::vector<PdbResidue>& segment : structure.segments) {
        std::vector<std::size_t> rows;
        for (std::size_t n = 0; n < segment.size(); ++n) {
            if (nucleotide(segment[n])) {
                rows.push_back(k + n);
            }
        }
        if (!rows.empty()) {
            out.push_back(std::move(rows));
        }
        k += segment.size();
    }
    return out;
}

std::vector<std::vector<std::size_t>> annotationStrandGroups(const std::vector<std::vector<std::size_t>>& strands,
                                                             const std::vector<ResiduePair>& pairs) {
    std::vector<std::size_t> parent(strands.size());
    for (std::size_t c = 0; c < parent.size(); ++c) {
        parent[c] = c;
    }
    std::map<std::int64_t, std::size_t> strandOf;
    for (std::size_t n = 0; n < strands.size(); ++n) {
        for (const std::size_t k : strands[n]) {
            strandOf[static_cast<std::int64_t>(k)] = n;
        }
    }
    const auto root = [&parent](std::size_t c) {
        while (parent[c] != c) {
            parent[c] = parent[parent[c]];
            c = parent[c];
        }
        return c;
    };
    for (const auto& [i, j] : pairs) {
        const auto si = strandOf.find(i);
        const auto sj = strandOf.find(j);
        if (si != strandOf.end() && sj != strandOf.end()) {
            const std::size_t a = root(si->second);
            const std::size_t b = root(sj->second);
            if (a != b) {
                parent[std::max(a, b)] = std::min(a, b);
            }
        }
    }
    std::vector<std::size_t> roots;  // dict insertion order
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t c = 0; c < strands.size(); ++c) {
        const std::size_t r = root(c);
        const auto at = std::find(roots.begin(), roots.end(), r);
        if (at == roots.end()) {
            roots.push_back(r);
            groups.push_back({c});
        } else {
            groups[static_cast<std::size_t>(at - roots.begin())].push_back(c);
        }
    }
    return groups;
}

AnnotationResult annotate(const PdbSource& source, const AnnotateOptions& options) {
    checkEngine(options);
    const std::string original = source.content();
    const PdbStructure s = io::readPdbText(original);
    Notes notes;
    for (const std::string& note : s.notes) {
        if (startsWith(note, "Kept model 1 of")) {
            notes.push_back(modelNote(note));
        } else if (startsWith(note, "Repaired")) {
            notes.emplace_back("The input has ChimeraX-damaged ATOM/HETATM lines; they are copied unchanged (repair "
                               "them with 'mdbridge pdb clean')");
        } else if (startsWith(note, "Dropped") || startsWith(note, "Ignored")) {
            if (startsWith(note, "Dropped the MASTER")) {
                notes.push_back(note);  // insertHeaderRecords drops it too
            }  // anything else is copied verbatim
        } else {
            notes.push_back(note);
        }
    }
    Annotated a = annotateStructure(s, std::move(notes), options);
    AnnotationResult result;
    result.text = io::insertHeaderRecords(original, a.records);
    result.notes = std::move(a.notes);
    result.records = std::move(a.records);
    result.engine = "builtin";
    return result;
}

AnnotationResult annotate(const PdbStructure& structure, const AnnotateOptions& options) {
    checkEngine(options);
    PdbStructure s = structure;
    Notes notes;
    for (const std::string& note : s.notes) {
        notes.push_back(startsWith(note, "Kept model 1 of") ? modelNote(note) : note);
    }
    Annotated a = annotateStructure(s, std::move(notes), options);
    std::erase_if(s.header, [](const std::string& line) { return startsWith(line, kAnnotationRemark); });
    s.header.insert(s.header.end(), a.records.begin(), a.records.end());
    AnnotationResult result;
    result.text = chem::writePdb(s, false);
    result.notes = std::move(a.notes);
    result.records = std::move(a.records);
    result.engine = "builtin";
    return result;
}

}  // namespace mdbridge::chem
