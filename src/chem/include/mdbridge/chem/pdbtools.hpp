#pragma once
// PDB hygiene for D- and L-RNA structures, as MDATK v0.0.4.1's PDB tools (mdatk/pdbtools.py;
// MDATK's own code, ported) on the fixed-column model of io/pdb.hpp: L-nucleotides from the signed
// chiral volumes of the sugar, residue ranges, the residue-name schemes (amber, ccd, rosetta, pdb),
// the ATOM / HETATM record policies, chain restore / split / combine / order, 5'-phosphate removal,
// renumbering, SEQRES and FASTA, the `pdb clean` chain of operations and `pdb check`'s summary.
//
// Every operation works in place on an io::PdbStructure and returns MDATK's notes (Notes, one line
// each, in MDATK's order and wording); impossible requests throw io::PdbError with MDATK's message.
// The outputs are byte for byte MDATK's (tests/chem/test_pdbtools.cpp, tests/data/chem/
// pdbtools_cases.json) except where MDBridge deliberately differs (docs/DIFFERENCES.md "Chemistry"):
//   CH-01 a residue number below -999 cannot be written: PdbError("residue number -1002 (A) does not
//         fit the PDB columns (-999 to 9999)") from writePdb (io::writePdb), nothing written
//   CH-03 combineChains(into) takes one character: PdbError("chain ids must be single characters")
//   CH-04 clean(stripConect = false) with repeated input serials says the CONECT records were dropped
//   CH-05 a residue with an RNA name on a 2'-deoxy sugar (C2', no O2') is DNA (MDATK 0.1.0)
//   CH-06 a chain end bonded O3'-P (<= 2.0 A) to a residue the schemes do not know (PSU, a 5'-GTP)
//         is no end (MDATK 0.1.0, without its LINK records)
//   CH-07 SEQRES names have three characters: longer nucleotide names as their CCD component, other
//         long names UNK, with notes (MDATK 0.0.17)
//   CH-13 PDB text is always text: a file is read through PdbSource::file (the command line names a
//         missing file); MDATK parsed a name that is not a file as PDB text
//   CH-15 a reversed residue range (A:4-1) is an error
//
// Callers: the `pdb clean / check / to-lrna` commands (milestone 6 W2-C), the annotation (W1-B:
// nucleotide, sequenceChains), the builder (W1-D: convertNames, setRecords, writePdb,
// sugarChirality, chiralVolume), the PDB clean-up window and milestone 7's generators (parseRanges,
// assignChirality, convertNames through chem/generator_api.hpp).

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "mdbridge/core/value.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

using Notes = std::vector<std::string>;
using Vec3 = std::array<double, 3>;
// atom name -> xyz in atom order, first atom of a name (io::PdbResidue::coords())
using AtomCoords = std::vector<std::pair<std::string, Vec3>>;

// SCHEMES, RECORD_POLICIES
inline constexpr std::array<std::string_view, 4> kSchemes = {"amber", "ccd", "rosetta", "pdb"};
inline constexpr std::array<std::string_view, 3> kRecordPolicies = {"keep", "rosetta", "pdb"};
// PHOSPHATE_5P: the atoms strip5primePhosphates removes ('*' read as ')
inline constexpr std::array<std::string_view, 7> kPhosphate5p = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"};
// MIN_CHIRAL_VOLUME (A^3): a centre below it is flat and prevents a call
inline constexpr double kMinChiralVolume = 0.5;
// CH-06: O3'-P distance (A) up to which a neighbour bonds to a chain end (MDATK 0.1.0)
inline constexpr double kO3PBonded = 2.0;

// --- residues ------------------------------------------------------------------------------------

// _nucleotide: io::parseResidueName of the name, guarded by a backbone atom (a primed name, '*'
// or P), so that a lanthanum ion LA is no nucleotide; with CH-05 an RNA name on a sugar with C2'
// and no O2' gives sugar "dna".
std::optional<io::ResidueName> nucleotide(const io::PdbResidue& residue);
// _is_polymer: a nucleotide or an amino acid (standard or Amber variant HID, CYX ...)
bool isPolymer(const io::PdbResidue& residue);
// _is_solvent: not a polymer residue, and water (HOH, WAT, TIP3 ...) or a single atom (an ion)
bool isSolvent(const io::PdbResidue& residue);
// _code: the one-letter code (base, amino acid) or 'X'
char residueCode(const io::PdbResidue& residue);
// _chirality: the assigned chirality, else the name's ("D" / "L"); nullopt for other residues
std::optional<std::string> residueChirality(const io::PdbResidue& residue);

// scheme_name: the name of an RNA nucleotide (base A C G U; chirality 'D' / 'L'; terminal "",
// "5", "3", "N") in `scheme`. Throws io::PdbError("unknown residue name scheme 'x' (choose from
// amber, ccd, rosetta, pdb)").
std::string schemeName(char base, char chirality, std::string_view terminal, std::string_view scheme);

// _terminals: the terminal state of every residue of a TER segment ("5", "3", "N" for a single
// nucleotide, "" inside and for residues that are no nucleotides), CH-06 included.
std::vector<std::string> segmentTerminals(const std::vector<io::PdbResidue>& segment);

// --- chirality -----------------------------------------------------------------------------------

// chiral_volume: a . (b x c) of the vectors from `centre` (CIF _chem_comp_chir sign); IEEE
// operations in MDATK's order.
double chiralVolume(const Vec3& centre, const Vec3& a, const Vec3& b, const Vec3& c);

// sugar_chirality: 'D', 'L' or no call from the volumes at C4', C3' and C1' (N9 for A / G, N1 for
// C / U / T; no base: N9 if present, else N1). No call when no centre is usable, centres disagree
// or one is flat. `volumes` in centre order, the usable centres only.
struct SugarChirality {
    std::optional<char> call;
    std::vector<std::pair<std::string, double>> volumes;
};
SugarChirality sugarChirality(const AtomCoords& coords, std::optional<char> base = std::nullopt);

// --- residue ranges ------------------------------------------------------------------------------

// (chain, start, end); start nullopt = the whole chain
struct ResidueRange {
    std::string chain;
    std::optional<std::int64_t> start;
    std::optional<std::int64_t> end;
    bool operator==(const ResidueRange&) const = default;
};

// What MDATK's functions accept as ranges: None, a text ("C:45-73, A:1-5", "C:45", "C", ":45-73";
// items separated by commas, semicolons or white space) or a list of texts and (chain, start, end)
// tuples (a tuple's end nullopt = its start; its chain stripped).
class RangeSpec {
public:
    using Item = std::variant<std::string, ResidueRange>;
    RangeSpec() = default;
    RangeSpec(std::nullopt_t) {}
    RangeSpec(std::string text) : m_value(std::move(text)) {}
    RangeSpec(const char* text) : m_value(std::string(text)) {}
    RangeSpec(std::vector<Item> items) : m_value(std::move(items)) {}

    bool isNone() const { return m_value.index() == 0; }
    // Python's truth value: false for None, "" and []
    bool truthy() const;
    const std::string* text() const { return std::get_if<std::string>(&m_value); }
    const std::vector<Item>* items() const { return std::get_if<std::vector<Item>>(&m_value); }

private:
    std::variant<std::monostate, std::string, std::vector<Item>> m_value;
};

// parse_ranges. Throws io::PdbError("cannot read residue range 'AB:1' (expected e.g. C:45-73, C:45
// or C)") and, CH-15, io::PdbError("residue range A:4-1 ends before it starts").
std::vector<ResidueRange> parseRanges(const RangeSpec& spec);
// _in_ranges
bool inRanges(const io::PdbResidue& residue, const std::vector<ResidueRange>& ranges);

// --- names and records ---------------------------------------------------------------------------

// assign_chirality: Residue::chirality of every nucleotide from the override ranges (forced L) >
// the sugar geometry (unless !detect) > the name; the notes of MDATK.
Notes assignChirality(io::PdbStructure& structure, bool detect = true, const RangeSpec& lRanges = {});

// convert_names: RNA nucleotides renamed to `scheme` (only those in `only` when given); terminals
// false: chain ends get the internal names. Uses the assigned chirality (else the name's).
Notes convertNames(io::PdbStructure& structure, std::string_view scheme, const RangeSpec& only = {},
                   bool terminals = true);

// set_records: "keep", "rosetta" (polymer residues ATOM) or "pdb" (names other than the standard
// ones HETATM). Throws io::PdbError("unknown record policy 'x' (choose from keep, rosetta, pdb)").
Notes setRecords(io::PdbStructure& structure, std::string_view policy);

// --- chains --------------------------------------------------------------------------------------

// restore_chains with a reference: chain ids (and with copyNumbers the residue numbers and
// insertion codes) copied residue by residue (water and ions excluded), matched by one-letter code.
Notes restoreChainsFromReference(io::PdbStructure& structure, const io::PdbStructure& reference,
                                 bool copyNumbers = false);
// restore_chains without one: `chainIds` per polymer TER segment (stripped, one character each),
// else A, B, C ...
Notes restoreChains(io::PdbStructure& structure,
                    const std::optional<std::vector<std::string>>& chainIds = std::nullopt);
// split_chains: `at` None: every further TER segment of a chain gets a new unused id; else the
// chains are split before the residues of `at` (file order; "C:45": residue 45 and the rest of C).
Notes splitChains(io::PdbStructure& structure, const RangeSpec& at = {});
// combine_chains (CH-03: `into` of more than one character is an error)
Notes combineChains(io::PdbStructure& structure, const std::vector<std::string>& chains,
                    const std::optional<std::string>& into = std::nullopt);
// order_chains: listed chains first in `order` (alphabetical without one), the others after them,
// water and ion segments last
Notes orderChains(io::PdbStructure& structure, const std::optional<std::vector<std::string>>& order = std::nullopt);

// --- numbering, phosphates, sequence -------------------------------------------------------------

// strip_5prime_phosphates: P / OP1 / OP2 / OP3 off the first nucleotide of every TER segment
// (lone phosphate residues removed; HO5' 0.96 A from O5' when the residue has hydrogens)
Notes strip5primePhosphates(io::PdbStructure& structure);
// renumber_residues
Notes renumberResidues(io::PdbStructure& structure, std::int64_t start = 1, bool perChain = false,
                       bool keepGaps = false);

// _sequence_chains: chain -> its residues in the polymer segments (water and ions excluded), in
// order of first appearance. The pointers point into `structure`.
std::vector<std::pair<std::string, std::vector<const io::PdbResidue*>>> sequenceChains(
    const io::PdbStructure& structure);
// CH-07 (MDATK 0.0.17 _seqres_names): the SEQRES names per chain, three characters each, and the
// notes that report the names it wrote in place of longer ones
struct SeqresNames {
    std::vector<io::SeqresChain> chains;
    Notes notes;
};
SeqresNames seqresNames(const io::PdbStructure& structure);
// seqres_records: the SEQRES lines of seqresNames (13 names per line)
std::vector<std::string> seqresRecords(const io::PdbStructure& structure);
// fasta: one entry per chain, its kind (RNA, DNA, L-RNA, protein ..., "L-nucleotides at 3-4")
std::string fasta(const io::PdbStructure& structure, std::string_view title = "structure");
// write_pdb(structure, strip_conect, seqres): io::writePdb with seqresNames' chains
std::string writePdb(const io::PdbStructure& structure, bool stripConect = true, bool seqres = false);

// summary (`pdb check`'s JSON, as a table in MDATK's key order): chains (residues, numbers,
// sequence, names, chirality counts), segments, other_residues, conect_records, notes
core::Value summary(const io::PdbStructure& structure);

// --- whole files ---------------------------------------------------------------------------------

// A PDB input: text, parsed as given, or a file (io::readPdbFile: std::filesystem::filesystem_error
// when it cannot be opened, core::ValueError for bytes that are not UTF-8). CH-13.
class PdbSource {
public:
    static PdbSource text(std::string text);
    static PdbSource file(std::filesystem::path path);
    io::PdbStructure read() const;
    // the PDB text: the text given, or the file's (io::readPdbFileText; the same errors as read)
    std::string content() const;
    // the FASTA title of clean without a title: the file's name without its extension; "structure"
    // for text
    std::string defaultTitle() const;

private:
    bool m_isFile = false;
    std::string m_text;
    std::filesystem::path m_path;
};

// clean's keyword arguments (MDATK's defaults)
struct CleanOptions {
    std::optional<std::string> scheme;   // None / "" = keep the names
    std::optional<std::string> records;  // None / "" = follows the scheme (rosetta, pdb, else keep)
    bool stripConect = true;
    bool detectL = true;
    RangeSpec lRanges;
    std::optional<PdbSource> chainReference;
    std::optional<std::vector<std::string>> chainIds;
    bool assignChains = false;
    bool copyNumbers = false;
    bool splitSegments = false;
    RangeSpec splitAt;
    std::vector<std::string> combine;
    std::optional<std::string> combineInto;
    // sort_chains: True = alphabetical; a non-empty chainOrder = that order (sort_chains=[...])
    bool sortChains = false;
    std::vector<std::string> chainOrder;
    bool renumber = false;
    std::int64_t start = 1;
    bool perChain = false;
    bool keepGaps = false;
    bool seqres = false;
    bool makeFasta = false;
    std::optional<std::string> title;
    bool terminalNames = true;
    bool strip5p = false;
};

struct CleanResult {
    std::string text;
    Notes notes;
    std::optional<std::string> fasta;
    io::PdbStructure structure;
};

// clean: read -> stale REMARK 950 dropped -> chain restore -> chirality -> split / combine / order
// -> 5'-phosphates -> names -> renumbering -> record policy -> CONECT / atom-count / SEQRES notes
// -> write (MDATK's order)
CleanResult clean(const PdbSource& source, const CleanOptions& options = {});

// to_lrna: L forced in `ranges`, only those renamed (the to_lrna.sh replacement)
CleanResult toLrna(const PdbSource& source, const RangeSpec& ranges, std::string_view scheme = "amber",
                   bool stripConect = true, bool terminals = true);

// `pdb check`: read, assignChirality(detect, lRanges), summary with the chirality notes after the
// read notes. The command prints it as json.dumps(info, indent=2) (core::toJsonAscii).
core::Value check(const PdbSource& source, const RangeSpec& lRanges = {});

// _listing: the first `limit` items joined with ", ", then " and N more"
std::string listing(const std::vector<std::string>& items, std::size_t limit = 8);

}  // namespace mdbridge::chem
