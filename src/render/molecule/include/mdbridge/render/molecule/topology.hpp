#pragma once
// Atom roles for the 3D representations (MDATK's mdatk/render/molecule/topology.py), computed once
// per figure from the cached `atoms` and `trace` analyses.
//
// Nucleotides are recognised by their atoms, not their residue names, so modified residues work
// when they keep the standard ring atoms:
//   * purine ring perimeter N9 C8 N7 C5 C6 N1 C2 N3 C4 (glycosidic N9),
//   * pyrimidine ring N1 C2 N3 C4 C5 C6 (glycosidic N1; C5 for pseudouridine),
//   * ribose C1' C2' C3' C4' O4'.
// Old atom names (C1*, O1P) are normalised. Amino acids use CA as their backbone point. Backbone
// segments break at chain changes, gaps in the residue numbering and control points more than
// 8.5 Å apart.

#include <array>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::render::molecule {

using Vec3 = std::array<double, 3>;

inline constexpr std::array<std::string_view, 9> kPurineRing = {"N9", "C8", "N7", "C5", "C6",
                                                                 "N1", "C2", "N3", "C4"};
inline constexpr std::array<std::string_view, 6> kPyrimidineRing = {"N1", "C2", "N3", "C4", "C5", "C6"};
inline constexpr std::array<std::string_view, 5> kSugarRing = {"C1'", "C2'", "C3'", "C4'", "O4'"};
inline constexpr double kSegmentBreak = 8.5;  // Å between consecutive control points

// "C1*" -> "C1'", "O1P" -> "OP1" (stripped of surrounding white space).
std::string normaliseName(std::string_view name);
// 'A', 'G', 'C', 'U' or 'T' for standard and common Amber / CHARMM residue names, else '?'.
char baseLetter(std::string_view resname);
// Bondi van der Waals radius (Å) of an upper-case element symbol; 1.7 for unknown elements.
double vdwRadius(std::string_view element);
// CPK-like colour (RGB 0..1) of an upper-case element symbol; (0.6, 0.4, 0.7) for unknown ones.
Vec3 elementRgb(std::string_view element);
// Colour of a nucleotide letter ('?' and unknown letters: grey).
Vec3 nucleotideRgb(char letter);

// What Topology is built from: the atoms analysis (per atom) and the trace analysis (per residue).
struct TopologyInput {
    std::vector<std::string> names;      // per atom
    std::vector<std::string> elements;   // per atom ("" -> carbon)
    std::vector<std::string> resnames;   // per atom
    std::vector<std::int32_t> row;       // per atom: its residue row (0 .. n_res-1)
    std::int32_t n_res = 0;
    std::vector<std::array<std::int32_t, 2>> bonds;  // atom index pairs
    std::vector<Vec3> ref_xyz;           // per atom: the reference frame (segment breaks)
    std::vector<std::int64_t> chain;     // per residue
    std::vector<std::int64_t> resindex;  // per residue (topology numbering)
    std::string backbone_atom = "P";     // "P" or "C4'" (anything else: "P")
};

// Per-residue atom indices (into the atoms arrays) for every style. -1 = none. The members are
// public and may be edited (as MDATK's tests do) before styles are built from a copy.
struct Topology {
    Topology() = default;
    explicit Topology(const TopologyInput& input);

    std::int32_t n_res = 0;
    std::vector<std::int32_t> row;
    std::vector<std::string> names;     // normalised
    std::vector<std::string> elements;  // stripped, upper case, "C" when blank
    std::vector<std::array<std::int32_t, 2>> bonds;
    std::vector<std::int64_t> chain;
    std::vector<std::int64_t> resindex;
    std::vector<std::map<std::string, std::int32_t, std::less<>>> table;  // per residue: name -> first atom
    std::vector<std::string> resname;  // per residue: its first atom's residue name, "?" without atoms
    std::vector<char> letter;          // per residue
    std::vector<double> radius;        // per atom
    std::vector<Vec3> element_rgb;     // per atom
    std::vector<bool> is_carbon, is_hydrogen;  // per atom
    std::vector<std::int32_t> control;      // spline control atom per residue
    std::vector<std::int32_t> anchor_atom;  // the control atom, else the residue's first atom
    std::vector<std::int32_t> side;         // toward the bases (C1', C4', CB or O)
    std::vector<std::int32_t> c1, c4;
    std::vector<std::int32_t> purine_res, purine_gly;
    std::vector<std::array<std::int32_t, 9>> purine_ring;
    std::vector<std::int32_t> pyrimidine_res, pyrimidine_gly;
    std::vector<std::array<std::int32_t, 6>> pyrimidine_ring;
    std::vector<std::int32_t> sugar_res;
    std::vector<std::array<std::int32_t, 5>> sugar_ring;
    std::vector<std::vector<std::int32_t>> segments;  // residue rows of each backbone segment (>= 2)

    // Residue rows that have a base ring (sorted, unique).
    std::vector<std::int32_t> nucleic() const;
    // First atom index of a residue's atom name, or -1.
    std::int32_t atomOf(std::int32_t residue, std::string_view name) const;

    // The segments for reference coordinates `xyz` (per atom).
    std::vector<std::vector<std::int32_t>> computeSegments(std::span<const Vec3> xyz) const;
};

}  // namespace mdbridge::render::molecule
