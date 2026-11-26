#pragma once
// Restraint dictionaries for the L-RNA residues of the `ccd` naming scheme, as MDATK's
// mdatk/dictionaries.py makes them (MDATK's own code, ported; docs/specs/m6-plan.md 1.4):
//
// * the components 0A 0C 0G 0U (the wwPDB CCD entries, from the generated reference tables of
//   reference_tables.hpp) and their eight 5' / 3' variants, with MDATK's isomeric L SMILES;
// * smilesMapping / checkSmiles / identifySmiles: a SMILES string against a component's atom-name
//   graph, hydrogens and sugar stereochemistry;
// * finishDictionary: an acedrg / eLBOW dictionary (mmCIF, io::parseCif) renamed to the CCD atom
//   names (heavy-atom graph isomorphism, hydrogens named from coordinates where two can swap),
//   chirality checked, leaving atoms removed, the component id / name / group set; written back
//   with io::writeCif, byte for byte MDATK's text.
//
// Errors are io::DictionaryError with MDATK's messages. Differences from MDATK: docs/DIFFERENCES.md
// "Chemistry: restraint dictionaries" (DIC-01 ...), CHD-01 (guanine's H21 / H22).

#include <array>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/chem/smiles.hpp"
#include "mdbridge/io/cif.hpp"

namespace mdbridge::chem {

// _chem_comp.group of a finished dictionary by default (Coot / Refmac / Phenix link RNA residues),
// and the groups the command line offers.
inline constexpr std::string_view kDefaultGroup = "RNA";
std::span<const std::string_view> dictionaryGroups();  // RNA, non-polymer

struct ComponentAtomName {
    std::string name;
    std::string element;
    int charge = 0;
};

struct ComponentBondRow {
    std::string atom1;
    std::string atom2;
    int order = 1;
};

// How the two hydrogens of a CH2 / NH2 group are told apart in coordinates: kind "+" / "-" is the
// sign of the chiral volume V(x; a, b, h), "cis" / "trans" the torsion h-x-a-b.
struct HydrogenRule {
    std::string hydrogen;
    std::string x;
    std::string a;
    std::string b;
    std::string kind;
};

// A sugar centre: sign of V(order[0], order[1], order[2]) about the centre; order[3] is its H.
struct ChiralityReference {
    std::string centre;
    std::array<std::string, 4> order;
    int sign = 1;
};

// One residue code: the molecule given to acedrg / eLBOW (atoms, bonds, hydrogens of `smiles`) and
// its reference atom names; `leaving` atoms are removed from the dictionary after renaming.
struct Component {
    std::string code;
    std::string parent;
    std::string base;      // A C G U
    std::string terminal;  // "" (CCD entry), "5" or "3"
    std::string name;      // _chem_comp.name
    std::string smiles;
    std::vector<ComponentAtomName> atoms;  // heavy atoms
    std::vector<ComponentBondRow> bonds;
    std::vector<std::pair<std::string, std::vector<std::string>>> hydrogens;  // heavy atom -> H names
    std::vector<HydrogenRule> hpos;
    std::vector<std::string> leaving;
    std::string description;

    // The hydrogen names of a heavy atom (empty for none).
    const std::vector<std::string>& hydrogenNames(std::string_view heavy) const;
    // Atom names of the finished dictionary's residue: heavy atoms first (CCD order), then hydrogens.
    std::vector<std::string> residueAtoms() const;
    std::string glycosidicN() const;  // N9 (purines) or N1
    std::vector<ChiralityReference> chirality() const;
    // The reference graph (names set); leaving = false: without the leaving atoms.
    MolGraph graph(bool withLeaving = true) const;
    std::string formula() const;
    bool isLeaving(std::string_view atom) const;
};

// COMPONENTS: 0A5 0A3 0C5 0C3 0G5 0G3 0U5 0U3 0A 0C 0G 0U.
const std::vector<Component>& dictionaryComponents();
// The component of exactly this code, nullptr if there is none.
const Component* findComponent(std::string_view code);
// get_component: the code stripped and upper-cased; DictionaryError("unknown residue code 'x'
// (known: 0A5, ...)").
const Component& getComponent(std::string_view code);

// describe_mismatch: "<what> (<formula>, N heavy atoms) does not match <code> (<formula>, M heavy
// atoms)".
std::string describeMismatch(const Component& comp, const MolGraph& tgt, std::string_view what);

struct SmilesMapping {
    std::map<int, std::string> names;  // SMILES heavy-atom index -> reference name
    std::vector<std::string> notes;
};

// smiles_mapping: the SMILES's atoms named after the component; refuses another compound, another
// protonation outside the phosphate and a missing or other configuration at C1' ... C4'.
SmilesMapping smilesMapping(const Component& comp, std::string_view smiles);
// check_smiles: the notes when the SMILES is a valid isomeric SMILES of the code.
std::vector<std::string> checkSmiles(std::string_view code, std::string_view smiles);
// identify_smiles: the codes whose generation molecule the SMILES encodes.
std::vector<std::string> identifySmiles(std::string_view smiles);

struct FinishOptions {
    const Component* reference = nullptr;  // default: the component of the code
    std::string group = std::string(kDefaultGroup);
    std::string name;      // _chem_comp.name (default: the component's)
    std::string sourceId;  // the component to take when the input holds several
    bool allowProtonation = false;
    bool keepNames = false;  // only set id / name / group
};

struct FinishResult {
    std::string text;
    std::vector<std::string> notes;
    // old atom name -> new name ("" = a removed leaving atom), in MDATK's order (heavy atoms in
    // the reference order, then hydrogens)
    std::vector<std::pair<std::string, std::string>> mapping;
    std::string sourceId;
    std::string program;  // "acedrg", "phenix.elbow" or ""
};

// finish_dictionary: the dictionary text renamed to `code` and its reference atom names.
FinishResult finishDictionary(std::string_view text, std::string_view code, const FinishOptions& options = {});
// rename_dictionary: the text of finishDictionary.
std::string renameDictionary(std::string_view text, std::string_view code, const FinishOptions& options = {});

}  // namespace mdbridge::chem
