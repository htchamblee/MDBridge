#pragma once
// Internal adapter: the two PDB-hygiene functions the helix builder needs (pdbtools.convert_names
// and set_records), until the PDB tools of milestone 6 W1-A (chem/pdbtools.hpp) land. At that merge
// the builder calls W1-A's functions instead and this file goes; it implements only what the builder
// reaches (no `only` ranges; every segment the builder makes is all nucleotides, so the CH-06
// neighbour test of _terminals never applies).

#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem::detail {

// convert_names(structure, scheme, only=None, terminals): renames the RNA residues, returns the notes.
std::vector<std::string> convertNames(io::PdbStructure& structure, std::string_view scheme, bool terminals);

// set_records(structure, policy) for "keep", "rosetta" and "pdb": returns the notes.
std::vector<std::string> setRecords(io::PdbStructure& structure, std::string_view policy);

// scheme_name(base, chirality, terminal, scheme)
std::string schemeName(char base, char chirality, std::string_view terminal, std::string_view scheme);

}  // namespace mdbridge::chem::detail
