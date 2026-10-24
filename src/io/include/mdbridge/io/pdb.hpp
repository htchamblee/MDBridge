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

// write_pdb(structure, strip_conect, seqres): "\n" line ends, every line right-stripped, final
// "\n". `seqres`: per-chain residue-name lists for SEQRES records (which residues: the caller's
// hygiene logic), inserted where MDATK inserts them.
struct SeqresChain {
    std::string chain;
    std::vector<std::string> names;
};
std::string writePdb(const PdbStructure& structure, bool stripConect = true,
                     const std::optional<std::vector<SeqresChain>>& seqres = std::nullopt);
std::vector<std::string> seqresLines(const std::vector<SeqresChain>& chains);

// pdbtools.parse_residue_name: nucleotide residue names (A, rG, LG5, DT, ADE ...).
struct ResidueName {
    char base = 'A';  // A C G T U
    std::string sugar;  // "rna" / "dna"
    char chirality = 'D';  // 'D' / 'L'
    std::string terminal;  // "", "5", "3", "N"
    bool operator==(const ResidueName&) const = default;
};
std::optional<ResidueName> parseResidueName(std::string_view name);

// annotate._insert_records: REMARK 950 and MASTER lines removed, `records` inserted before the
// first MODEL / ATOM / HETATM line, everything else verbatim.
std::string insertHeaderRecords(std::string_view text, const std::vector<std::string>& records);

// builder.transform_pdb_text: x' = R x + t on ATOM / HETATM lines (ANISOU rotated, half-even
// rounding), all else verbatim with its line ends; std::nullopt for unparsable coordinates.
// Throws PdbError("the fitted coordinates do not fit the PDB columns").
std::optional<std::string> transformPdbText(std::string_view text, const std::array<std::array<double, 3>, 3>& rotation,
                                            const std::array<double, 3>& translation);

// builder._tail: columns 28-76 of a new atom. Throws PdbError("coordinates do not fit the PDB
// columns").
std::string newAtomTail(const std::array<double, 3>& xyz, double occupancy = 1.0, double bfactor = 0.0);

}  // namespace mdbridge::io
