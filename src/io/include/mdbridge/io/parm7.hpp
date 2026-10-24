#pragma once
// Amber parm7 / prmtop topologies: MDATK's own section reader and writer (mdatk/io/parm7.py), the
// Lennard-Jones table its energy analysis reads, and the topology MDAnalysis' TOPParser builds
// (which is what MDATK's analyses and exports see).

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "mdbridge/core/topology.hpp"

namespace mdbridge::io {

// The values of one %FLAG section: strings (Fortran a), integers (i) or reals (e, f).
using Parm7Values = std::variant<std::vector<std::string>, std::vector<std::int64_t>, std::vector<double>>;
using Parm7Sections = std::vector<std::pair<std::string, Parm7Values>>;  // file order; a repeated flag: the last

// read_sections(path, wanted): every (or every wanted) section, decoded with its %FORMAT. Throws
// core::ValueError with MDATK's / Python's messages ("<path>: cannot parse format line '...'",
// "invalid literal for int() with base 10: '...'", "could not convert string to float: '...'"),
// std::filesystem::filesystem_error when the file cannot be read.
Parm7Sections readParm7Sections(const std::filesystem::path& path,
                                const std::optional<std::set<std::string>>& wanted = std::nullopt);

// The section `flag`, or nullptr.
const Parm7Values* findSection(const Parm7Sections& sections, std::string_view flag);

// Lennard-Jones coefficients of atom pairs (parm7 LENNARD_JONES_ACOEF / BCOEF).
struct LjTable {
    std::vector<std::int64_t> typeIndex;  // 1-based type of every atom
    std::int64_t ntypes = 0;
    std::vector<std::int64_t> nbIndex;
    std::vector<double> acoef;
    std::vector<double> bcoef;

    // Throws core::ValueError("<path>: missing parm7 sections ['NONBONDED_PARM_INDEX', ...]").
    static LjTable fromParm7(const std::filesystem::path& path);
    // (A, B) of each pair (i[k], j[k]) of 0-based atoms; 0 for pairs without a 6-12 term.
    // Throws std::out_of_range for indices outside the tables.
    std::pair<std::vector<double>, std::vector<double>> pairCoefficients(std::span<const std::int64_t> i,
                                                                        std::span<const std::int64_t> j) const;
};

// One section to write: kind 'a' (20a4), 'i' (10I8), 'e' (5E16.8) or 's' (1a80).
struct Parm7Section {
    std::string flag;
    char kind = 'i';
    Parm7Values values;  // strings for a / s, integers for i, reals for e
};

// _format_values: the lines of one section.
std::vector<std::string> formatParm7Values(char kind, const Parm7Values& values);

// write_parm7: the fixed %VERSION line, the TITLE section, then the sections; LF line ends on every
// platform. Throws core::ValueError for non-ASCII text, std::filesystem::filesystem_error.
void writeParm7(const std::filesystem::path& path, const std::vector<Parm7Section>& sections,
                std::string_view title = "MDATK");

// The topology MDAnalysis' TOPParser builds from a parm7 file: names, types, charges (float32
// arithmetic / 18.2223), masses, elements from ATOMIC_NUMBER, resnames, resids 1..n, residues
// from RESIDUE_POINTER, one segment SYSTEM or RESIDUE_CHAINID segments with per-atom chain IDs,
// ids 1..n, bonds of both bond sections. No coordinates. Throws core::ValueError with MDAnalysis'
// messages ("<path> is not a valid TOP file. %VE Missing in header", ...).
std::shared_ptr<core::Topology> readParm7Topology(const std::filesystem::path& path);

// .parm7 / .prmtop / .top, case-insensitive (os.path.splitext rules).
bool isParm7(std::string_view path);

}  // namespace mdbridge::io
