#pragma once
// MDATK's structure builder (mdatk/builder.py v0.0.4.1, MDATK's own code, ported), part 3: base
// mutations. A base is swapped for another keeping the glycosidic frame (a weighted Kabsch fit of the
// new base's template on C1' and the glycosidic N and its two ring neighbours), the sugar-phosphate
// atoms verbatim; with `compensatory` the Watson-Crick / G·U partner of every mutated base (found by
// the builtin base-pair detector in the input) gets the complement too. The PDB text is MDATK's byte
// for byte (the Kabsch rotation comes from MDBridge's own 3 x 3 SVD, so a third decimal may differ
// from NumPy / LAPACK's in rare cases: the parity tests then compare the coordinates within 0.0011 A),
// with MDBridge's fixes (docs/DIFFERENCES.md):
//   CH-08 a mutation without a chain id reads the whole number: "15G" = blank chain, residue 15
//         (MDATK 0.1.0; MDATK 0.0.4.1 read chain "1", residue 5)
//   CH-09 (D11) a mutated base gets ideal hydrogens (the builder templates' hydrogens, placed by the
//         same fit as its heavy atoms, after them in the hydrogen-rule order) when the old base had
//         any: the note "B:3: replaced 4 hydrogen(s) of the old base (H8 H61 H62 H2) with 3 of the
//         new base (H3 H5 H6) at ideal positions" (MDATK: the old base's hydrogens dropped, none added)
//
// Callers: `build mutate` (milestone 6 W2-C), the builder window (W3-B), the planner's mutation step
// (milestone 7). Errors are io::BuildError with MDATK's messages (io::PdbError from reading and
// writing the PDB text).

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

// One base change (builder.Mutation)
struct Mutation {
    std::string chain;  // "" = blank chain
    std::int64_t resseq = 0;
    std::string icode;  // "" or one letter
    char base = 'A';    // A C G U
    bool isExplicit = true;  // false: a compensatory change of the partner
    std::string reason;      // "compensates A:5 G" for a compensatory change

    // "<chain or _>:<resseq><icode>"
    std::string where() const;
    bool operator==(const Mutation&) const = default;
};

// parse_mutations: "A:5G, B:20C" (items separated by commas, semicolons or white space; also
// "A:5>G", "A:5->G", "A:5:G", "A:5AG" with insertion code A, ":5G" and, CH-08, "5G" for the blank
// chain; the base case-insensitive), or a list of such items. BuildError "cannot read mutation 'x'
// (expected e.g. A:5G - chain, residue number, new base)", "residue A:5 is mutated twice", "no
// mutations given (e.g. A:5G, B:20C)".
std::vector<Mutation> parseMutations(std::string_view spec);
std::vector<Mutation> parseMutations(std::span<const std::string> items);

// mutated_name: the residue name for `base` in the scheme, chirality and terminal state of `name`
// (LG5 -> LC5, 0A3 -> 0G3, GUA -> URA). BuildError "cannot rename residue 'x' (RNA nucleotides only)".
std::string mutatedName(std::string_view name, char base);

struct MutateOptions {
    // also change the Watson-Crick / G·U partner (detected in the input) of every mutated base to
    // the complement, unless the partner is mutated explicitly
    bool compensatory = true;
    // with `compensatory`: a partner that forms a G·U pair with the new base is kept
    bool wobble = false;
    double cutoff = 3.5;  // the base-pair detector's H-bond cutoff (A)
};

struct MutateResult {
    std::string text;  // PDB text (chem::writePdb: CONECT stripped, REMARK 950 dropped)
    Notes notes;       // the read notes, the compensation notes, one or more per changed base, ...
};

// mutate(source, mutations, compensatory, wobble, cutoff): swap the bases; the notes in MDATK's
// order ("A:2 G->U (G -> U, glycosidic-frame fit 0.12 Å)", "... - compensatory, compensates A:2 U",
// "Dropped 2 REMARK 950 base-pair annotation lines ...", "Stripped 1 CONECT records", "Mutated 2
// residue(s)", then the chirality WARNINGs). BuildError "residue A:99 not found", "B:2 DG is not an
// RNA nucleotide", "A:2 G: missing atom(s) C8 of the glycosidic frame", "...: the base swap changed
// the sugar chirality call (D -> L, {...})"; the mutation-spec errors of parseMutations.
MutateResult mutate(const PdbSource& source, std::string_view mutations, const MutateOptions& options = {});
MutateResult mutate(const PdbSource& source, std::vector<Mutation> mutations, const MutateOptions& options = {});

}  // namespace mdbridge::chem
