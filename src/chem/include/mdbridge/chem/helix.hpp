#pragma once
// MDATK's structure builder (mdatk/builder.py), part 1: the Arnott A-RNA templates with their ideal
// hydrogens, the sequence functions and the levelled messages the builder window and the planner
// show (check_helix), and the Watson-Crick helix (build_helix with pairing "wc", no triplex): a
// duplex or a single strand of either hand, byte for byte as MDATK writes it, with MDBridge's fixes
// (docs/DIFFERENCES.md CH-01, CH-02, CHD-01, D13 "REMARK 951 MDBridge IDEAL").
//
// Python's None is std::nullopt where MDATK tells it apart from "" (partner, pairing, triplex,
// third, scheme). Errors are io::BuildError (io::PdbError for the CH-01 column check), with MDATK's
// messages.

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

using io::BuildError;

// Arnott A-RNA screw: 11 residues per turn
inline constexpr double kHelixTwist = 32.7;  // degrees per residue
inline constexpr double kHelixRise = 2.81;   // Angstrom per residue

// The pairings (PAIRINGS) and triplex motifs (TRIPLEX_MOTIFS) MDATK knows; only "wc" without a
// triplex is built here (the others: chem/helix_fit, milestone 6 W2-A).
inline constexpr std::array<std::string_view, 5> kPairings = {"wc", "reverse-wc", "hoogsteen", "reverse-hoogsteen",
                                                              "hh"};
inline constexpr std::array<std::string_view, 2> kTriplexMotifs = {"hoogsteen", "reverse-hoogsteen"};
inline constexpr std::array<std::string_view, 4> kSchemes = {"amber", "ccd", "rosetta", "pdb"};

// --- templates -----------------------------------------------------------------------------------

using Vec3 = std::array<double, 3>;

struct TemplateAtom {
    std::string name;
    Vec3 xyz;
};

// template(base, hydrogens): the A-form monomer of base A, C, G or U (case ignored) in the helix
// frame (D sugar, strand 1, level 0): the 12 backbone heavy atoms, the base heavy atoms (CCD order),
// then with `hydrogens` the sugar and base hydrogens (BOND_H lengths; guanine's H21 trans and H22 cis
// to N1 as the CCD names them, CHD-01). Throws BuildError("no template for base 'X' (A, C, G, U)").
std::vector<TemplateAtom> helixTemplate(std::string_view base, bool hydrogens = true);

// BASE_NAMES[base] / the hydrogen names of _H_RULES[base] (empty for another letter) and of the sugar.
std::span<const std::string_view> baseAtomNames(char base);
std::span<const std::string_view> baseHydrogenNames(char base);
std::span<const std::string_view> sugarHydrogenNames();

// --- sequences -----------------------------------------------------------------------------------

struct ParsedSequence {
    std::string sequence;
    std::vector<std::string> notes;
};
// parse_sequence: whitespace, 5'- / -3' and '-' ignored, T read as U (with a note).
ParsedSequence parseSequence(std::string_view text, std::string_view what = "sequence");
// complement: the Watson-Crick complement 5'->3' of a cleaned sequence.
std::string complement(std::string_view sequence);

// check_pairing (None / "" -> "wc") and check_triplex (None / "" / "none" / "no" -> nullopt).
std::string checkPairing(std::optional<std::string_view> pairing);
std::optional<std::string> checkTriplex(std::optional<std::string_view> triplex);
// helix_kind: "A-form", "parallel tWW", "hoogsteen triplex" ...
std::string helixKind(std::optional<std::string_view> pairing = "wc", std::optional<std::string_view> triplex = {});
// auto_partner of a cleaned sequence for a checked pairing.
std::string autoPartner(std::string_view sequence, std::string_view pairing = "wc");

struct ResolvedPartner {
    std::string sequence;
    std::optional<std::string> partner;  // 5'->3'; nullopt: a single strand
    std::vector<std::string> notes;
};
// resolve_partner(sequence, partner, pairing): partner nullopt or "auto", "none", or a sequence.
ResolvedPartner resolvePartner(std::string_view sequence, std::optional<std::string_view> partner = "auto",
                               std::optional<std::string_view> pairing = "wc");

struct ResolvedThird {
    std::string sequence;
    std::vector<std::string> notes;
};
// resolve_third(sequence, third, motif) of a cleaned sequence and a checked motif.
ResolvedThird resolveThird(std::string_view sequence, std::optional<std::string_view> third = "auto",
                           std::string_view motif = "hoogsteen");

// pair_types: "wc", "wobble" or "mismatch" for (seq[i], partner[n-1-i]).
std::vector<std::string> pairTypes(std::string_view sequence, std::string_view partner);

struct PairKind {
    char first;    // strand-1 base
    char partner;  // partner (or third-strand) base
    std::string kind;  // wc / wobble / mismatch (wc), canonical / mismatch (others)
    bool operator==(const PairKind&) const = default;
};
std::vector<PairKind> pairingKinds(std::string_view sequence, std::string_view partner,
                                   std::optional<std::string_view> pairing = "wc");
std::vector<PairKind> triplexKinds(std::string_view sequence, std::string_view third, std::string_view motif);

// One canonical H-bond of a pair (_PAIR_HBONDS): the atom of the strand-1 (duplex) base, the atom of
// the partner, which of them donates (1 or 2), the donor's hydrogen ("H3+": N3 proton of C+). With
// CHD-01 the guanine N2 donor of wc and reverse-wc G·C and of reverse-Hoogsteen G·G is H22.
struct PairHbond {
    std::string_view atom1;
    std::string_view atom2;
    int donor;
    std::string_view hydrogen;
};
// The H-bonds of `rule` (a pairing or triplex motif) between bases a and b; empty: not a canonical pair.
std::span<const PairHbond> pairHbonds(std::string_view rule, char a, char b);

struct HelixMessage {
    std::string level;  // "error", "warning", "info"
    std::string text;
    bool operator==(const HelixMessage&) const = default;
};
// check_helix: the sequence-level messages of a build request (every pairing and triplex).
std::vector<HelixMessage> checkHelix(std::string_view sequence, std::optional<std::string_view> partner = "auto",
                                     std::optional<std::string_view> pairing = "wc",
                                     std::optional<std::string_view> triplex = {},
                                     std::optional<std::string_view> third = "auto");

// _parse_ids(value, count, str / int): "A,B" -> {A, B}, "A" -> {A, A}; a triplex's missing third
// value is the next free chain id (C, D ...) / residue 1. Throws BuildError.
std::vector<std::string> parseChainIds(std::string_view value, int count);
std::vector<std::int64_t> parseStartNumbers(std::string_view value, int count);

// --- Watson-Crick helices ------------------------------------------------------------------------

struct HelixOptions {
    std::optional<std::string> partner = "auto";  // "auto", "none" or the partner 5'->3' (nullopt: auto)
    std::string chirality = "D";                  // D / L (case and spaces ignored; "" -> D)
    std::string chains = "A,B";                   // chain ids of strand 1, 2 (one value: both)
    std::string start = "1,1";                    // first residue numbers of strand 1, 2
    std::optional<std::string> scheme;            // amber / ccd / rosetta / pdb; nullopt or "": A C G U / 0A ...
    bool terminals = true;                        // 5' / 3' names in the schemes that have them
    bool hydrogens = false;
    bool phosphate5 = false;                      // P, OP1, OP2, OP3 on each 5' end (default 5'-OH)
    bool center = true;                           // helix centre at the origin (axis along z)
    std::optional<std::string> pairing;           // nullopt / "" / "wc" (others: helix_fit, W2-A)
    std::optional<std::string> triplex;           // nullopt / "" / "none" (triplexes: helix_fit, W2-A)
    std::optional<std::string> third = "auto";    // the third strand of a triplex
};

// build_helix for Watson-Crick duplexes and single strands: the PDB text; the build notes are
// appended to `notes` (when given) only on success. A pairing other than wc or a triplex throws
// BuildError("... is built by the fitted helix builder (not available yet)") until helix_fit lands.
std::string buildHelix(std::string_view sequence, const HelixOptions& options = {},
                       std::vector<std::string>* notes = nullptr);

}  // namespace mdbridge::chem
