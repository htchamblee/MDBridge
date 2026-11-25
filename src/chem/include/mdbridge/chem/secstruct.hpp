#pragma once
// Nucleic-acid secondary structure as MDATK computes it (mdatk/secstruct.py): dot-bracket notation
// and the builtin base-pair detection. Ported from MDATK's own code; results are bit for bit
// MDATK's (the same IEEE operations in numpy's order, tests/chem/test_secstruct.cpp).
//
// Dot-bracket: dotbracketToPairs reads (), [], {}, <> and Aa ... Zz pseudoknot levels with '&' or
// whitespace (Python's str.isspace) between chains; pairsToDotbracket writes it, putting each helix
// on the lowest bracket level it does not cross (longest helices first).
//
// Builtin detection (PairDetector): canonical pairs from the donor-acceptor distances of their
// base-pair hydrogen bonds (pairRules()), each <= cutoff, with the C1' atoms 8.5-12.5 A apart
// (residues without a C1' / C1* atom skip that test). Every base pairs at most once; candidates
// are taken greedily by the smallest sum of the two shortest H-bond distances rounded to the
// 0.05 A grid, near-ties going to the C1'-C1' distance closer to 10.5 A.
//
// Not ported: MDATK's BarnabaRunner (the GPL engine barnaba is never used by MDBridge).

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mdbridge::chem {

// Pair kinds (the values of base-pair maps): WATSON_CRICK, WOBBLE, OTHER.
inline constexpr int kWatsonCrick = 1;
inline constexpr int kWobble = 2;
inline constexpr int kOther = 3;
// KIND_LABELS: "Watson-Crick", "wobble", "non-canonical"; "" for any other value.
std::string_view kindLabel(int kind);

inline constexpr std::string_view kOpeners = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
inline constexpr std::string_view kClosers = ")]}>abcdefghijklmnopqrstuvwxyz";
inline constexpr std::string_view kUnpaired = ".-,:_";
inline constexpr char kChainSeparator = '&';

// PAIR_RULES: (base 1, partner bases, kind, (atom of base 1, atom of partner) ..., bonds required).
struct PairRule {
    char base;
    std::string_view partners;
    int kind;
    std::vector<std::pair<std::string_view, std::string_view>> atoms;
    int required;
};
// G-C (3 bonds, 2 required), A-U / A-T (2 of 2), G-U / G-T wobble (2 of 2), in that order.
const std::vector<PairRule>& pairRules();

inline constexpr double kC1WindowLow = 8.5;  // C1_WINDOW (A)
inline constexpr double kC1WindowHigh = 12.5;
inline constexpr double kC1Ideal = 10.5;   // C1_IDEAL
inline constexpr double kScoreTie = 0.05;  // SCORE_TIE: the grid of the H-bond score (A)
// BASE_PAIR_ATOMS: every atom name of pairRules(), sorted.
const std::vector<std::string>& basePairAtoms();

// base_letter: the base (A C G U T) of a nucleotide residue name in any scheme MDATK reads
// (io::parseResidueName), else nullopt.
std::optional<char> baseLetter(std::string_view resname);
// canonical_name: the standard PDB name (A C G U; DA DC DG DT), else nullopt.
std::optional<std::string> canonicalName(std::string_view resname);
// lw_kind: pair kind of a Leontis-Westhof class ("WCc" Watson-Crick, "GUc" wobble, else other).
int lwKind(std::string_view code);

// --- dot-bracket

using ResiduePair = std::pair<std::int64_t, std::int64_t>;

struct DotBracket {
    std::vector<ResiduePair> pairs;  // (i, j), i < j, 0-based residue indices, sorted
    std::int64_t n = 0;              // residues (separators do not count)
    std::vector<std::int64_t> chainBreaks;  // first residue of every chain after the first
    bool operator==(const DotBracket&) const = default;
};

// dotbracket_to_pairs of UTF-8 text (characters are code points, as Python counts them; a byte
// that is not valid UTF-8 is a character of its own, shown as \udcXX). Throws core::ValueError
// with MDATK's messages ("unbalanced dot-bracket: ...", "invalid dot-bracket character ...").
DotBracket dotbracketToPairs(std::string_view text);

// _stems: pairs grouped into helices of directly stacked pairs (i, j), (i + 1, j - 1) ...,
// starting from the pairs sorted.
std::vector<std::vector<ResiduePair>> stems(std::span<const ResiduePair> pairs);
// _crosses: a and b form a pseudoknot (a0 < b0 < a1 < b1 or b0 < a0 < b1 < a1).
bool crosses(const ResiduePair& a, const ResiduePair& b);

// pairs_to_dotbracket: dot-bracket text of n residues, '&' before every index of chainBreaks
// inside (0, n). Pairs are normalised (i < j) and deduplicated. Throws core::ValueError for a
// pair outside the residues or with i == j ("invalid base pair (a, b) for n residues"), a residue
// in two different pairs, and more than 30 bracket levels.
std::string pairsToDotbracket(std::span<const ResiduePair> pairs, std::int64_t n,
                              std::span<const std::int64_t> chainBreaks = {});

// unique_pairs: pairs where every residue occurs once (first come) and the rejected rest.
struct UniquePairs {
    std::vector<ResiduePair> keep;
    std::vector<ResiduePair> dropped;
};
UniquePairs uniquePairs(std::span<const ResiduePair> pairs);

// --- builtin detection

// One row of PairDetector.candidates: coordinate-array residues `first` (base of the rule) and
// `second` (partner), the rule's kind, the H-bond score (sum of the two shortest donor-acceptor
// distances) and |C1'-C1' - 10.5| (10.5 when a residue has no C1' atom).
struct PairCandidate {
    std::int64_t first = 0;
    std::int64_t second = 0;
    int kind = 0;
    double score = 0.0;
    double c1Deviation = 0.0;
};

// One detected pair: residues i < j and the kind.
struct BasePair {
    std::int64_t i = 0;
    std::int64_t j = 0;
    int kind = 0;
    bool operator==(const BasePair&) const = default;
};

// Atom name -> row of the coordinate array, per residue.
using AtomIndex = std::map<std::string, std::int64_t, std::less<>>;

class PairDetector {
public:
    // resnames: residue name per residue; atomIndex: per residue, atom name -> row of the
    // coordinates later passed to candidates / detect (negative rows count as missing). Residues
    // that are no nucleotide, or lack the atoms of a rule, never pair by that rule.
    PairDetector(std::span<const std::string> resnames, std::span<const AtomIndex> atomIndex,
                 double cutoff = 3.5, bool wobble = true);

    // Maps an atom-atom difference vector (first minus second) to the vector whose length is the
    // distance, e.g. its periodic minimum image.
    using DifferenceImage = std::function<std::array<double, 3>(const std::array<double, 3>&)>;

    // Every pair meeting a rule, in MDATK's order: rule by rule, first residues in order, partners
    // in order. xyz: 3 values per row (x y z), float64. Throws std::out_of_range when a row the
    // rules read is not in xyz. With `image`, every distance the rules read (donor-acceptor and
    // C1'-C1') is the length of image(p - q) instead of |p - q| (MDBridge's periodic base pairs,
    // docs/DIFFERENCES.md AN-11).
    std::vector<PairCandidate> candidates(std::span<const double> xyz, const DifferenceImage* image = nullptr) const;
    // The pairs, sorted by (i, j): candidates taken greedily (see select).
    std::vector<BasePair> detect(std::span<const double> xyz, const DifferenceImage* image = nullptr) const;

    // detect's greedy step on given candidates: stable order by (round(score / 0.05) * 0.05 with
    // round half to even, then c1Deviation; NaN last), a candidate is kept when neither residue is
    // used yet; result sorted. Residue indices must be in [0, nResidues).
    static std::vector<BasePair> select(std::span<const PairCandidate> candidates, std::size_t nResidues);

    std::size_t residueCount() const { return m_n; }
    double cutoff() const { return m_cutoff; }

private:
    struct Rule {
        int kind = 0;
        int required = 0;
        std::size_t width = 0;  // atoms per residue (bonds of the rule)
        std::vector<std::int64_t> rowsP;  // residues of base 1 that can reach `required` bonds
        std::vector<std::int64_t> atomsP;  // [rowsP.size() * width], -1 for missing
        std::vector<std::int64_t> c1P;
        std::vector<std::int64_t> rowsQ;
        std::vector<std::int64_t> atomsQ;
        std::vector<std::int64_t> c1Q;
    };
    std::size_t m_n = 0;
    double m_cutoff = 3.5;
    std::vector<Rule> m_rules;
};

}  // namespace mdbridge::chem
