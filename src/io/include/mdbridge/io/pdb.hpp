#pragma once
// PDB files as MDATK's PDB tools read and write them (mdatk/pdbtools.py read_pdb / write_pdb): a
// fixed-column model that keeps columns 28-76 of every atom record verbatim, repairs ChimeraX-damaged
// lines, keeps the first model, and writes TER after every segment; plus the small PDB text
// utilities of the annotation and builder code. Byte-exact with MDATK.

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::io {

// pdbtools.PDBError
class PdbError : public core::ValueError {
public:
    using core::ValueError::ValueError;
    const char* pythonName() const noexcept override { return "PDBError"; }
};

// builder.BuildError (a PDBError)
class BuildError : public PdbError {
public:
    using PdbError::PdbError;
    const char* pythonName() const noexcept override { return "BuildError"; }
};

// 62 chain IDs: A-Z, a-z, 0-9
inline constexpr std::string_view kChainIds = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

struct PdbAtom {
    std::string record;  // "ATOM" / "HETATM" as read (case kept)
    std::string name;
    std::string alt;
    std::string tail;  // columns 28-76 as read (49 characters)
    std::string element;  // columns 77-78 stripped
    std::string charge;  // columns 79-80 stripped
    std::optional<std::int64_t> serial;  // input serial, for CONECT remapping
    std::string anisou;  // columns 28-80 of the following ANISOU record, padded

    // float(tail[3:11]), float(tail[11:19]), float(tail[19:27]). Throws core::ValueError.
    std::array<double, 3> xyz() const;
    bool operator==(const PdbAtom&) const = default;
};

struct PdbResidue {
    std::string name;
    std::string chain;
    std::int64_t resseq = 0;
    std::string icode;
    std::vector<PdbAtom> atoms;
    std::optional<std::string> chirality;

    // "<chain or _>:<resseq><icode> <name>"
    std::string label() const;
    // name (with * as ') -> xyz of the first atom of that name, in atom order
    std::vector<std::pair<std::string, std::array<double, 3>>> coords() const;
    bool operator==(const PdbResidue&) const = default;
};

struct PdbStructure {
    std::vector<std::string> header;
    std::vector<std::vector<PdbResidue>> segments;
    std::vector<std::string> trailer;
    std::vector<std::string> conect;
    std::string model;  // the first MODEL line, "" without models
    std::vector<std::string> notes;

    std::vector<const PdbResidue*> residues() const;
    std::vector<std::string> chainIds() const;  // in order of first appearance
    bool operator==(const PdbStructure&) const = default;
};

// read_pdb of PDB text, with every note and PdbError message of MDATK.
PdbStructure readPdbText(std::string_view text);
// read_pdb of a file: strict UTF-8, universal newlines. Throws std::filesystem::filesystem_error,
// core::ValueError for undecodable bytes ("'utf-8' codec can't decode byte 0xe9 in position 10:
// invalid continuation byte").
PdbStructure readPdbFile(const std::filesystem::path& path);
// The text of a PDB file as readPdbFile reads it (strict UTF-8, the same errors), line ends as
// stored (every consumer splits lines as Python's splitlines, for which they are alike).
std::string readPdbFileText(const std::filesystem::path& path);

// write_pdb(structure, strip_conect, seqres): "\n" line ends, every line right-stripped, final
// "\n". `seqres`: per-chain residue-name lists for SEQRES records (which residues and names: the
// caller's hygiene logic, chem::seqresNames), inserted where MDATK inserts them. CH-01: a residue
// number below -999 does not fit columns 23-26: PdbError("residue number -1002 (A) does not fit the
// PDB columns (-999 to 9999)"), MDATK wrote "A-1002" and shifted the line; numbers above 9999 wrap
// as MDATK's do. CH-07: SEQRES names have at most three characters (PdbError otherwise; MDATK
// wrote longer names in four columns).
struct SeqresChain {
    std::string chain;
    std::vector<std::string> names;
};
std::string writePdb(const PdbStructure& structure, bool stripConect = true,
                     const std::optional<std::vector<SeqresChain>>& seqres = std::nullopt);
std::vector<std::string> seqresLines(const std::vector<SeqresChain>& chains);

// pdbtools.parse_residue_name: nucleotide residue names (A, rG, LG5, DT, ADE ...) and, as MDATK
// 0.1.x, L-DNA: the CCD's L-deoxyribonucleotides 0DA 0DC 0DG 0DT (data/sources/ccd) and MDATK's
// own LDA LDC LDG LDT (LDA5 / LDA3 / LDAN at chain ends in older files). A name alone can be
// something else (in the CCD, LA and RA are the lanthanum and radium ions, LDA a detergent): readers
// that see the atoms take a name for a nucleotide only with hasBackboneAtoms.
struct ResidueName {
    char base = 'A';  // A C G T U
    std::string sugar;  // "rna" / "dna"
    char chirality = 'D';  // 'D' / 'L'
    std::string terminal;  // "", "5", "3", "N"
    // CH-22: the CCD code of a modified ribonucleotide read as its parent (CH, the N3-protonated
    // cytidine: base C, D-RNA; core's kIndModifiedRibonucleotides), "" for the others
    std::string component;
    bool operator==(const ResidueName&) const = default;
};
std::optional<ResidueName> parseResidueName(std::string_view name);

// has_backbone_atoms (MDATK 0.1.x): whether a residue with these atom names has a sugar-phosphate
// atom (a primed name, an old-style '*' or P), the guard of parseResidueName's names.
bool isBackboneAtomName(std::string_view atomName);

// Python's str.upper() (full case mappings: "ß" -> "SS"), as the PDB tools compare residue names.
std::string pythonUpper(std::string_view text);
template <class Names>
bool hasBackboneAtoms(const Names& atomNames) {
    for (const auto& name : atomNames) {
        if (isBackboneAtomName(name)) {
            return true;
        }
    }
    return false;
}

// standard_atom_name (MDATK 0.1.x): the standard (PDB v3) name of atom `name` of residue `resname`.
// The CCD names some atoms of 0DA and 0DT its own way (0DA's N9A ... HN61, H2P; 0DT's C5M:
// data/sources/ccd); those give the standard name (N9, H61, HOP2, C7), every other atom its own.
std::string standardAtomName(std::string_view resname, std::string_view name);

// annotate._insert_records: REMARK 950 and MASTER lines removed, `records` inserted before the
// first MODEL / ATOM / HETATM line, everything else verbatim.
std::string insertHeaderRecords(std::string_view text, const std::vector<std::string>& records);

// builder.transform_pdb_text: x' = R x + t on ATOM / HETATM lines (ANISOU rotated, half-even
// rounding), all else verbatim with its line ends; std::nullopt for unparsable coordinates.
// Throws BuildError("the fitted coordinates do not fit the PDB columns").
std::optional<std::string> transformPdbText(std::string_view text, const std::array<std::array<double, 3>, 3>& rotation,
                                            const std::array<double, 3>& translation);

// builder._tail: columns 28-76 of a new atom. Throws BuildError("coordinates do not fit the PDB
// columns").
std::string newAtomTail(const std::array<double, 3>& xyz, double occupancy = 1.0, double bfactor = 0.0);

}  // namespace mdbridge::io
