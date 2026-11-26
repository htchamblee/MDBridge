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

// Arnott monomer: backbone heavy atoms (P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1').
std::span<const CylindricalAtom> arnottBackbone();
// Arnott monomer: the heavy atoms of base A, C, G or U (empty for any other letter), CCD order.
std::span<const CylindricalBaseAtom> arnottBase(char base);

// The 4JRD nucleotide (heavy atoms, CCD order) and its screw: one level on is a turn by
// -jrdTwist() degrees about z and a move of -jrdRise() Angstrom along z.
std::span<const CylindricalAtom> jrdNucleotide();
double jrdTwist();
double jrdRise();

std::span<const ComponentInfo> components();
// The rows of one component (code 0A ...), in the standard residue's order; empty for an unknown code.
std::span<const ComponentAtom> componentAtoms(std::string_view code);
std::span<const ComponentBond> componentBonds(std::string_view code);
std::span<const ComponentHydrogen> componentHydrogens(std::string_view code);
// The rules of one residue (A C G U or 0A 0C 0G 0U).
std::span<const HydrogenPosition> hydrogenPositions(std::string_view residue);
std::span<const ChiralCentre> lRiboseChirality();

}  // namespace mdbridge::chem::reference
