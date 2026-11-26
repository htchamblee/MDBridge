#pragma once
// Base-pair and dot-bracket annotation of a PDB file as REMARK 950 records, as MDATK v0.0.4.1's
// mdatk/annotate.py writes it (MDATK's own code, ported). annotate reads a PDB input, finds the base
// pairs of its first model with the builtin canonical-pair detector (chem::PairDetector, the
// `basepairs` analysis' detector) and inserts the records before the first MODEL / ATOM / HETATM
// record; every other line (all models, atom serials, CONECT ...) is copied unchanged, the REMARK
// 950 lines of an earlier annotation are replaced and a MASTER record is dropped
// (io::insertHeaderRecords). Record layout (fixed columns; a residue field is the chain id, the
// residue number in 4 columns, the insertion code and the name left-justified in 3 columns):
//
//     REMARK 950 MDBridge BASE-PAIR ANNOTATION
//     REMARK 950 ENGINE builtin, canonical pairs from H-bond distances <= 3.5 A
//     REMARK 950 PAIRS 20 (WATSON-CRICK 18, WOBBLE 2, OTHER 0)
//     REMARK 950 PAIR  A  12 G   B  33 C   WCc          one per pair (WCc / GUc)
//     REMARK 950 CHAINS A,B                             one block per group of strands joined by pairs
//     REMARK 950 SEQ   GGGCAAAGCU&AGCUUUGCCC            '&' between strands
//     REMARK 950 DBN   ((((...(((&)))...))))            Watson-Crick and wobble pairs
//
// SEQ and DBN are wrapped to 60 characters per line (continuation lines repeat the keyword). A strand
// is the nucleotides (chem::nucleotide) of one TER segment, so files with blank chain ids still get
// a '&' between strands. Milestone 7's refinement reads the PAIR records by keyword.
//
// Where MDBridge deliberately differs (docs/DIFFERENCES.md "Chemistry"; tests/chem/
// test_annotate.cpp, tests/data/chem/annotate_cases.json):
//   D13   the first record names MDBridge (MDATK: "REMARK 950 MDATK BASE-PAIR ANNOTATION")
//   CH-14 no Barnaba (D9): engine "auto" is the builtin detector; "barnaba" runs the builtin detector
//         with the note "WARNING: Barnaba is not available in MDBridge; using the builtin detector"
//         before the engine's notes; the stacking note has no installation hint; no STACK records
//   CH-13 a file is read through PdbSource::file; text is never taken for a path
// Messages that name the program say 'mdbridge pdb clean'.

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

// ENGINES
inline constexpr std::array<std::string_view, 3> kAnnotationEngines = {"auto", "builtin", "barnaba"};
inline constexpr std::string_view kAnnotationRemark = "REMARK 950";
inline constexpr std::size_t kAnnotationWrap = 60;  // characters of SEQ / DBN per record
// CH-14
inline constexpr std::string_view kBarnabaNote =
    "WARNING: Barnaba is not available in MDBridge; using the builtin detector";

// annotate's keyword arguments (MDATK's defaults). cutoff / wobble: the builtin detector's (as the
// basepairs analysis); stacking only decides the "no stacking" note (the builtin engine finds none).
struct AnnotateOptions {
    bool pairs = true;
    bool stacking = true;
    bool dotbracket = true;
    std::string engine = "auto";
    double cutoff = 3.5;
    bool wobble = true;
};

struct AnnotationResult {
    std::string text;                  // the annotated PDB file
    Notes notes;
    std::vector<std::string> records;  // the REMARK 950 lines written
    std::string engine;                // the engine used: always "builtin"
};

// annotate of PDB text or a file: the records inserted into the input's text (io::insertHeaderRecords).
// Throws core::ValueError("unknown annotation engine 'x' (choose from auto, builtin, barnaba)")
// before reading, io::PdbError from the reader, core::ValueError for unreadable coordinates, and
// for a file what PdbSource::read throws.
AnnotationResult annotate(const PdbSource& source, const AnnotateOptions& options = {});
// annotate of a Structure: its header without REMARK 950 lines plus the records, written by
// writePdb(stripConect = false); the structure's notes are reported as they are.
AnnotationResult annotate(const io::PdbStructure& structure, const AnnotateOptions& options = {});

// _res_field: chain (' ' when blank, first character), number right-justified in 4, insertion code
// (' ' when blank), name left-justified in 3
std::string annotationResidueField(const io::PdbResidue& residue);
// _wrapped: "REMARK 950 KEYWD text" per 60 characters (one record for empty text)
std::vector<std::string> annotationWrapped(std::string_view keyword, std::string_view text);
// _strands: the indices (into structure.residues()) of the nucleotides of every TER segment that
// has any
std::vector<std::vector<std::size_t>> annotationStrands(const io::PdbStructure& structure);
// _strand_groups: strand numbers grouped by the pairs between them (file order within and between
// groups)
std::vector<std::vector<std::size_t>> annotationStrandGroups(const std::vector<std::vector<std::size_t>>& strands,
                                                             const std::vector<ResiduePair>& pairs);

}  // namespace mdbridge::chem
