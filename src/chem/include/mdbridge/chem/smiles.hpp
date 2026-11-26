#pragma once
// Molecular graphs, SMILES and graph matching as MDATK's restraint dictionaries use them (the
// "molecular graphs", "SMILES" and "graph matching" parts of mdatk/dictionaries.py; MDATK's own
// code, ported). Errors are io::DictionaryError with MDATK's messages; text is UTF-8 and counted
// in characters (code points), as Python counts a str.

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/io/cif.hpp"

namespace mdbridge::chem {

using io::DictionaryError;

// MolGraph: heavy atoms with hydrogen counts; bond order 0 = aromatic / delocalised / unknown.
struct MolGraph {
    // A neighbour in a stereo centre's SMILES order: an atom index, or kHydrogen for a hydrogen
    // (the bracket H or an explicit [H] atom folded into the centre).
    static constexpr int kHydrogen = -1;
    struct Stereo {
        std::vector<int> order;  // all four neighbours, in SMILES order
        std::string tag;         // "@" or "@@"
    };

    std::vector<std::string> elements;
    std::vector<std::optional<int>> charges;  // nullopt: unknown (a partial charge in a dictionary)
    std::vector<int> hcount;
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> hnames;  // hydrogen names per heavy atom (dictionaries)
    // neighbours and bond orders, in the order the bonds were added (Python's dict order: a bond
    // added again keeps its place and takes the new order)
    std::vector<std::vector<std::pair<int, int>>> adj;
    std::map<int, Stereo> stereo;  // SMILES tetrahedral centres

    int addAtom(std::string element, std::optional<int> charge = 0, int hydrogens = 0, std::string name = {},
                std::vector<std::string> hydrogenNames = {});
    void addBond(int a, int b, int order = 1);
    // Removes the bond a-b from both atoms' lists (no-op without it).
    void removeBond(int a, int b);
    std::size_t size() const { return elements.size(); }
    // The order of the bond a-b, nullopt without one.
    std::optional<int> bond(int a, int b) const;
    // Hill formula including hydrogens, with the net charge when every charge is known and the
    // sum is not zero: "C10H12N5O8P-2", "ClH4N".
    std::string formula() const;
};

// parse_smiles: heavy-atom graph of a SMILES string (the first whitespace-separated word), with
// hydrogen counts (implicit ones of the organic subset from the lowest fitting valence, aromatic
// atoms one less), charges and tetrahedral stereo tags (@, @@, @TH1, @TH2). Explicit [H] atoms
// with one neighbour are folded into it. Throws DictionaryError ("empty SMILES", "SMILES: cannot
// read 'x' at position N", "SMILES: unclosed ring bond 1, 2", ...).
MolGraph parseSmiles(std::string_view smiles);

// At most this many isomorphisms are enumerated (MDATK's MAX_MATCHES).
inline constexpr std::size_t kMaxMatches = 2000;

// _all_matches: every heavy-atom isomorphism ref -> tgt (element + connectivity) in MDATK's
// backtracking order, up to kMaxMatches; mapping[i] is the tgt index of ref atom i.
std::vector<std::vector<int>> allMatches(const MolGraph& ref, const MolGraph& tgt);

// _score: hydrogen-count (4), charge (2) and bond-order (1) disagreements of a mapping; unknown
// charges and orders (0) are not counted.
int matchScore(const MolGraph& ref, const MolGraph& tgt, const std::vector<int>& mapping);

// match_graph: the mapping of lowest score (ties: the lexicographically smallest mapping),
// nullopt when the graphs are not isomorphic.
std::optional<std::vector<int>> matchGraph(const MolGraph& ref, const MolGraph& tgt);

}  // namespace mdbridge::chem
