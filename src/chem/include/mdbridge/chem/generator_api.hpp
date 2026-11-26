#pragma once
// The chem functions milestone 7's Lua script generators call while they check and write plans
// (docs/specs/m6-plan.md 4.4), as value-in / value-out functions in one header that milestone 7
// binds to Lua without wrapping logic: strings, lists and small records in, the same out; notes as
// string lists; errors as io::PdbError / io::BuildError / io::DictionaryError (and
// std::filesystem::filesystem_error for files) with MDATK's messages.
//
// What MDATK v0.0.4.1's generators use, and where it is here:
//   workflows/amber.py      parse_ranges, read_pdb, assign_chirality        parseRanges, readModel
//   workflows/buildstep.py  resolve_partner, resolve_third, pair_types,     resolvePartner, resolveThird,
//                           check_helix, helix_kind, read_pdb,               pairTypes, checkHelix, helixKind,
//                           parse_mutations                                  readModel; parse_mutations
//                                                                           with milestone 6 W2-B (chem/mutate)
//   workflows/dictstep.py   read_pdb, assign_chirality, convert_names       readModel(scheme "ccd"),
//                           (ccd), COMPONENTS, DEFAULT_GROUP,                components, kDefaultGroup,
//                           check_smiles, parse_smiles                       checkSmiles, parseSmiles;
//                           fetch_pubchem_smiles                             the program's, not chem's (D7:
//                                                                           app/commands/pubchem.hpp)
//   workflows/refine.py     read_pdb, _nucleotide, assign_chirality,         readModel (nucleotide per
//                           read_annotation (REMARK 950 PAIR / STACK),       residue), readAnnotation,
//                           secstruct.lw_kind                                lwKind
//
// The embedded copies of pdbtools.py / dictionaries.py in generated scripts are milestone 7's
// decision (plan 4.5); the `mdbridge pdb clean` / `dict finish` commands keep MDATK's options so
// that either answer works.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/chem/dictionaries.hpp"
#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/chem/smiles.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem::generator {

// ---- PDB models ---------------------------------------------------------------------------------

// parse_ranges: "C:45-73, A:1-5", "C", ":45-73" (nullopt / "": none). io::PdbError.
std::vector<ResidueRange> parseRanges(const std::optional<std::string>& text);

// One residue of a model as the generators look at it.
struct ModelResidue {
    std::string chain;
    std::int64_t number = 0;
    std::string icode;
    std::string name;   // as read, or after the scheme of readModel
    std::string label;  // "A:5 G"
    // _nucleotide: base, sugar, the name's chirality and terminal (nullopt for other residues)
    std::optional<io::ResidueName> nucleotide;
    // the assigned chirality "D" / "L" (readModel with `assign`); else the name's for
    // nucleotides; nullopt for other residues
    std::optional<std::string> chirality;
    std::vector<std::string> atoms;  // atom names in file order
};

struct Model {
    std::vector<ModelResidue> residues;  // file order
    std::vector<std::string> header;     // the header records (REMARK 950 ...)
    Notes readNotes;                     // read_pdb's notes (ChimeraX repair, models ...)
    Notes chiralityNotes;                // assign_chirality's notes (WARNING: ... lines included)
    Notes nameNotes;                     // convert_names' notes (with a scheme)
};

// read_pdb(path), then assign_chirality(detect=True, l_ranges) when `assign`, then
// convert_names(scheme) when a scheme is given. Throws filesystem_error for a missing file,
// io::PdbError for a file that is no PDB and for bad ranges or schemes.
Model readModel(const std::filesystem::path& path, bool assign = true,
                const std::optional<std::string>& lRanges = std::nullopt,
                const std::optional<std::string>& scheme = std::nullopt);

// The dictionary codes (COMPONENTS order; only the 5' / 3' variants unless !variantsOnly) that name
// residues of `model` (dictstep.py: the restraint dictionaries a model needs, with scheme "ccd").
std::vector<std::string> codesInModel(const Model& model, bool variantsOnly = true);

// ---- REMARK 950 (annotate.hpp writes them) ------------------------------------------------------

struct RemarkResidue {
    std::string chain;
    std::int64_t number = 0;
    std::string icode;
    std::string name;
    bool operator==(const RemarkResidue&) const = default;
};

struct RemarkPair {
    RemarkResidue a;
    RemarkResidue b;
    std::string code;  // WCc, GUc ... (lwKind)
    bool operator==(const RemarkPair&) const = default;
};

struct Annotation {
    bool annotated = false;          // any REMARK 950 record
    std::vector<RemarkPair> pairs;   // REMARK 950 PAIR
    std::vector<RemarkPair> stacks;  // REMARK 950 STACK (MDBridge writes none, CH-14)
};

// read_annotation(header) of MDATK's refinement step: the PAIR and STACK records (residue fields
// at columns 18-26 and 28-36, the code after column 37); unreadable records are skipped.
Annotation readAnnotation(std::span<const std::string> header);

// secstruct.lw_kind
using chem::lwKind;

// ---- the builder's sequence functions (chem/helix.hpp) --------------------------------------------

using chem::checkHelix;
using chem::complement;
using chem::helixKind;
using chem::HelixMessage;
using chem::pairTypes;
using chem::parseSequence;
using chem::resolvePartner;
using chem::ResolvedPartner;
using chem::resolveThird;
using chem::ResolvedThird;

// ---- restraint dictionaries (chem/dictionaries.hpp) -------------------------------------------------

using chem::checkSmiles;
using chem::Component;
using chem::kDefaultGroup;
using chem::parseSmiles;

// COMPONENTS in MDATK's order
inline const std::vector<Component>& components() { return dictionaryComponents(); }
// Component.is_variant: a 5' / 3' variant (0A5 ...), not a CCD parent (0A ...)
inline bool isVariant(const Component& component) { return !component.terminal.empty(); }

}  // namespace mdbridge::chem::generator
