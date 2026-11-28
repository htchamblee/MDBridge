#pragma once
// Reference data of the structure builder and the restraint dictionaries, generated from public
// sources by tools/gen_chem_tables.py (src/chem/helix_tables.inc, component_tables.inc; provenance
// in data/sources/ATTRIBUTION.md, docs/DIFFERENCES.md CHD-01 / CHD-02):
//
// * the Arnott A-RNA fibre monomer (Arnott, Hukins & Dover 1972; AmberTools' nucgen.dat $ARNA) in
//   cylindrical coordinates of the helix frame (D sugar, strand 1, level 0), as MDATK's builder.py
//   `_BACKBONE` / `_BASES`;
// * the 4JRD nucleotide (PDB 4JRD chain A residue 10) in the frame of the screw of that duplex, with
//   the screw's twist and rise (MDATK: `_JRD_NUCLEOTIDE`, `JRD_TWIST`, `JRD_RISE`; CHD-02);
// * the experimental geometry of the Hoogsteen duplex and the pyrimidine-motif triplex (CH-18): the
//   mean U·A-U and C+·G-C base triples of the RNA triple helices of PDB 3P22, 4PLX, 6SVS, 7JNH in the
//   helix frame of their screw, the screws of both kinds with the spreads of the experimental steps,
//   the pairs' H-bond distances and the partner nucleosides' pose spreads (src/chem/triplex_tables.inc);
// * the wwPDB CCD's L-ribonucleotides 0A 0C 0G 0U as MDATK's dictionaries.py `COMPONENTS` describes
//   the parents (standard atom names, phosphate deprotonated), how the two hydrogens of a CH2 / NH2
//   group are named from coordinates (CHD-01) and the signs of the L-ribose chiral volumes.

#include <span>
#include <string_view>

namespace mdbridge::chem::reference {

struct CylindricalAtom {
    std::string_view name;
    double r;    // Angstrom from the helix axis
    double phi;  // degrees
    double z;    // Angstrom along the axis
};

struct CylindricalBaseAtom {
    char base;  // A C G U
    std::string_view name;
    double r;
    double phi;
    double z;
};

// One atom of a mean experimental base triple: level "AUU" (U·A-U) or "GCC" (C+·G-C), strand 0
// purine, 1 Watson-Crick partner, 2 Hoogsteen (third-strand) partner; Cartesian, Angstrom.
struct TriplexAtom {
    std::string_view level;
    int strand;
    std::string_view name;
    double x;
    double y;
    double z;
};

// The experimental screw of a kind ("triplex-hoogsteen", "hoogsteen"): origin and unit axis (one level
// on is a turn of -twist degrees about the axis and a move of -rise Angstrom along it), the standard
// deviations of the experimental steps' axis offsets (2, Angstrom), tilts (2, degrees), twist and rise
// about it (tools/gen_chem_tables.py step_parameters), the mean separation of consecutive base planes
// and its standard deviation, the number of steps and of independent steps (sites: the copies of one
// molecule in a crystal counted once; the standard deviations are over the sites), the inclination of the
// levels' base planes to the steps' axes and its standard deviation.
struct TriplexScrew {
    std::string_view kind;
    double ox, oy, oz;
    double ax, ay, az;
    double twist;
    double rise;
    double sd[6];
    double stack;
    double stackSd;
    int steps;
    int sites;
    double inclination;    // degrees: the levels' base planes against the steps' axes
    double inclinationSd;
    // D19: the extremes of the experimental steps (each step measured on its levels' C1' and common base
    // atoms): twist (degrees), rise (Angstrom), inclination (degrees), and per strand (the kind's strand
    // order; 0 past its strands) the largest stacking of consecutive bases (Angstrom)
    double twistRange[2];
    double riseRange[2];
    double inclinationRange[2];
    double stackMax[3];
};

// Donor-acceptor distance of a U·A-U pair ("wc" / "hoogsteen" rule; purine atom, partner atom): mean and
// standard deviation over the independent triples, their number.
struct TriplexHbond {
    std::string_view rule;
    std::string_view atom1;
    std::string_view atom2;
    double mean;
    double sd;
    int sites;
};

// RMS per component of a U·A-U partner's rotation (degrees) and C1' shift (Angstrom) about its mean
// pose relative to the purine nucleoside, by table strand (1, 2).
struct TriplexPoseSpread {
    int strand;
    double rotation;
    double shift;
    int sites;  // independent triples
};

struct ComponentInfo {
    std::string_view code;  // 0A 0C 0G 0U
    char base;
    std::string_view name;  // _chem_comp.name
};

struct ComponentAtom {
    std::string_view code;
    std::string_view name;
    std::string_view element;
    int charge;
};

struct ComponentBond {
    std::string_view code;
    std::string_view atom1;
    std::string_view atom2;
    int order;  // 1 single, 2 double
};

struct ComponentHydrogen {
    std::string_view code;
    std::string_view heavy;
    std::string_view hydrogen;
};

// kind "+" / "-": sign of the chiral volume V(x; a, b, h); "cis" / "trans": torsion h-x-a-b.
struct HydrogenPosition {
    std::string_view residue;  // A C G U (D) or 0A 0C 0G 0U (L)
    std::string_view hydrogen;
    std::string_view x;
    std::string_view a;
    std::string_view b;
    std::string_view kind;
};

// Sign of V(centre; a, b, c) of the L-ribose; "N" stands for N9 (purines) or N1 (pyrimidines).
struct ChiralCentre {
    std::string_view centre;
    std::string_view a;
    std::string_view b;
    std::string_view c;
    int sign;
};

// A CCD component that is a standard base with one proton more (the N3-protonated cytidine CH, the C+ of
// the builder's G·C+ Hoogsteen pairs): its code, the base, the extra hydrogen and its heavy atom.
struct ProtonatedBase {
    std::string_view code;
    char base;
    std::string_view proton;
    std::string_view heavy;
};

// Arnott monomer: backbone heavy atoms (P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1').
std::span<const CylindricalAtom> arnottBackbone();
// Arnott monomer: the heavy atoms of base A, C, G or U (empty for any other letter), CCD order.
std::span<const CylindricalBaseAtom> arnottBase(char base);

// The 4JRD nucleotide (heavy atoms, CCD order) and its screw: one level on is a turn by
// -jrdTwist() degrees about z and a move of -jrdRise() Angstrom along z.
std::span<const CylindricalAtom> jrdNucleotide();
double jrdTwist();
double jrdRise();

// The mean triple atoms of a level ("AUU" / "GCC") and strand (0 .. 2), in the tables' order.
std::span<const TriplexAtom> triplexLevel(std::string_view level, int strand);
// The screw of "triplex-hoogsteen" or "hoogsteen" (nullptr for another kind).
const TriplexScrew* triplexScrew(std::string_view kind);
std::span<const TriplexHbond> triplexHbonds();
std::span<const TriplexPoseSpread> triplexPoseSpreads();

std::span<const ComponentInfo> components();
// The rows of one component (code 0A ...), in the standard residue's order; empty for an unknown code.
std::span<const ComponentAtom> componentAtoms(std::string_view code);
std::span<const ComponentBond> componentBonds(std::string_view code);
std::span<const ComponentHydrogen> componentHydrogens(std::string_view code);
// The rules of one residue (A C G U or 0A 0C 0G 0U).
std::span<const HydrogenPosition> hydrogenPositions(std::string_view residue);
std::span<const ChiralCentre> lRiboseChirality();
// The CCD's N3-protonated cytidine (CH, HN3 on N3; tools/gen_chem_tables.py protonated_cytidine).
const ProtonatedBase& protonatedCytidine();

}  // namespace mdbridge::chem::reference
