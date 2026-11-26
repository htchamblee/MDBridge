#pragma once
// Residue groups from the chains of a topology (MDATK 0.0.10 / 0.1.1 mdatk/chaingroups.py, ported;
// owner decision D15 of milestone 5). No Qt.
//
// One group per polymer chain, water and ions left out. The chains are, as for the "chain" colour
// mode: the chain IDs when they name several chains (PDB, or a parm7 with RESIDUE_CHAINID), else the
// segment ids when they name several, else the bonded fragments in topology order (from the
// topology's bonds, or from consecutive residues within covalent distance when it has none). Every
// group's residue spec is checked to select exactly its chain's atoms with the selection code the
// analyses use: plain residue ranges when the residue numbers are unique, else chain-qualified
// ranges ("A:1-40"), else a selection.
//
// mergeChainGroups() decides what may change in a project. Chain groups only fill groups that are
// empty or, when a topology was just chosen (`automatic`), replace the untouched default groups
// (Aptamer / Target) when there are at least as many chains and the defaults do not already select
// a chain each; then every reference to a default follows the chain group that takes its place
// (core::renameGroup). Otherwise existing groups, the defaults included, are never changed: on
// request (the page's *Groups from chains*) the chains that no group selects yet are added.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::app {

// A residue group for one chain, checked against the topology it was made from.
struct ChainGroup {
    core::ResidueGroup group;
    std::string label;      // chain ID, segment id or fragment number ("1", "2", ...)
    std::string source;     // "chainID", "segid" or "fragment"
    std::size_t residues = 0;
    std::string coveredBy;  // an existing group that already selects exactly this chain's atoms
};

// What mergeChainGroups() changed.
struct ChainMerge {
    std::vector<std::string> added;                          // chain groups added
    std::vector<std::pair<std::string, std::string>> replaced;  // default group -> chain group in its place
    std::vector<std::string> references;                     // fields rewritten for the replaced defaults
    std::vector<std::string> covered;                        // chains an existing group already selects
    std::vector<std::string> pending;                        // chains left to the button (user-made groups)

    bool changed() const { return !added.empty() || !replaced.empty(); }
    // One line for the window; `source`: where the chains come from (a file name).
    std::string summary(const std::string& source) const;
};

// Per residue: neither water nor an ion (water by name or composition: one heavy atom, an oxygen).
std::vector<bool> polymerMask(const core::Topology& topology);

// (label, source, residue indices) of the polymer chains in topology order. `positions`: the
// coordinates (3 per atom) for linking residues of a topology without bonds; empty when it has none.
struct PolymerChain {
    std::string label;
    std::string source;
    std::vector<std::int32_t> residues;
};
std::vector<PolymerChain> polymerChains(const std::shared_ptr<const core::Topology>& topology,
                                        const std::vector<float>& positions);

// The groups are still exactly the defaults a new project starts with.
bool untouchedDefaults(const std::vector<core::ResidueGroup>& groups);

// One verified group per polymer chain. Names (Chain_A, or Chain_1 ... for fragments) and colours
// differ from those of `existing`; coveredBy names the existing group that already selects exactly a
// chain's atoms. `current`: the coordinates the selections may use (null for none).
std::vector<ChainGroup> chainGroups(const std::shared_ptr<const core::Topology>& topology, const core::Frame* current,
                                    const std::vector<core::ResidueGroup>& existing = {});

// chainGroups() of a topology or structure file read on its own (a parm7 has no coordinates).
// Throws what the readers throw.
std::vector<ChainGroup> topologyChainGroups(const std::filesystem::path& path,
                                            const std::vector<core::ResidueGroup>& existing = {});

// Brings the chain groups into `project` (see the top of this file).
ChainMerge mergeChainGroups(core::Project& project, const std::vector<ChainGroup>& chains, bool automatic = false);

}  // namespace mdbridge::app
