#pragma once
// The molecular system MDATK works on, as MDAnalysis models it: atoms -> residues -> segments,
// per-format optional attributes, atom groups, frames and trajectories.
//
// MDATK's behaviour depends on which attributes a format provides (chain IDs, elements, charges,
// bonds, insertion codes), so absent attributes are std::nullopt rather than defaults, and the io
// readers fill the model with MDAnalysis' conventions (see the io headers). Positions and boxes are
// float32, as MDAnalysis keeps them; MDATK widens them to double itself.

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::core {

// Distinct strings of a topology, so that string matching runs once per distinct value.
class StringPool {
public:
    StringPool();  // id 0 is ""
    std::uint32_t intern(std::string_view text);
    std::optional<std::uint32_t> find(std::string_view text) const;
    const std::string& operator[](std::uint32_t id) const { return m_strings[id]; }
    std::size_t size() const { return m_strings.size(); }

private:
    std::vector<std::string> m_strings;
    std::map<std::string, std::uint32_t, std::less<>> m_ids;
};

using Bond = std::array<std::int32_t, 2>;  // atom indices, first < second

// Atoms, residues and segments in file order, 0-based. Fill the vectors, then call finish().
// String attributes are ids into `strings`.
struct Topology {
    StringPool strings;

    // ---- atoms (one entry per atom)
    std::vector<std::uint32_t> names;
    std::vector<std::uint32_t> types;  // "" when the format has none
    std::vector<std::int32_t> atomResidue;  // residue index of every atom
    // the topology value as MDAnalysis stores it ("Mg" from PDB columns; "" for unknown entries)
    std::optional<std::vector<std::uint32_t>> elements;
    std::optional<std::vector<std::uint32_t>> chainIDs;  // per atom (PDB / mmCIF; parm7 RESIDUE_CHAINID)
    std::optional<std::vector<double>> masses;
    std::optional<std::vector<double>> charges;  // e (parm7 CHARGE / 18.2223)
    std::optional<std::vector<std::int64_t>> ids;  // atom serials of the file (the "id" keyword)
    std::optional<std::vector<std::uint32_t>> altLocs;
    std::optional<std::vector<std::uint32_t>> recordTypes;  // "ATOM" / "HETATM"
    std::optional<std::vector<double>> occupancies;
    std::optional<std::vector<double>> tempFactors;
    std::optional<std::vector<std::int64_t>> formalCharges;

    // ---- residues (one entry per residue)
    std::vector<std::int64_t> resids;
    std::vector<std::int64_t> resnums;
    std::vector<std::uint32_t> resnames;
    std::optional<std::vector<std::uint32_t>> icodes;  // PDB / mmCIF, even when all ""
    std::vector<std::int32_t> residueSegment;

    // ---- segments
    std::vector<std::uint32_t> segids;

    // canonical (first < second), sorted, unique; std::nullopt: the format has no bond information
    std::optional<std::vector<Bond>> bonds;

    // Checks that the sizes agree and indices are in range (std::invalid_argument otherwise), sorts
    // and deduplicates the bonds, and builds the residue -> atoms map.
    void finish();

    std::size_t atomCount() const { return names.size(); }
    std::size_t residueCount() const { return resids.size(); }
    std::size_t segmentCount() const { return segids.size(); }

    // Atoms of residue r in ascending index order (after finish()).
    std::span<const std::int32_t> residueAtoms(std::int32_t residue) const;

    const std::string& name(std::int32_t atom) const { return strings[names[atom]]; }
    const std::string& resname(std::int32_t residue) const { return strings[resnames[residue]]; }
    const std::string& segid(std::int32_t segment) const { return strings[segids[segment]]; }
    std::int32_t segmentOfAtom(std::int32_t atom) const { return residueSegment[atomResidue[atom]]; }

private:
    std::vector<std::int32_t> m_residueAtomOffsets;  // CSR, residueCount() + 1 entries
    std::vector<std::int32_t> m_residueAtoms;
};

// Connected components of the bond graph, numbered by increasing lowest atom index (MDAnalysis
// fragindices). Throws NoDataError when the topology has no bonds.
std::vector<std::int32_t> fragmentIndices(const Topology& topology);

// A list of atoms of one topology. Order and duplicates are significant (an MDAnalysis AtomGroup);
// selections always give sorted unique atoms.
class AtomGroup {
public:
    AtomGroup() = default;
    AtomGroup(std::shared_ptr<const Topology> topology, std::vector<std::int32_t> indices);
    // All atoms.
    explicit AtomGroup(std::shared_ptr<const Topology> topology);

    const Topology& topology() const { return *m_topology; }
    const std::shared_ptr<const Topology>& sharedTopology() const { return m_topology; }
    const std::vector<std::int32_t>& indices() const { return m_indices; }
    std::size_t size() const { return m_indices.size(); }
    bool empty() const { return m_indices.empty(); }
    std::int32_t operator[](std::size_t i) const { return m_indices[i]; }

    // ag[mask]: the atoms whose mask entry is true, order kept.
    AtomGroup masked(const std::vector<bool>& mask) const;
    // ag[positions]: positions into this group (negative counts from the end), order and
    // duplicates kept. Throws std::out_of_range.
    AtomGroup at(std::span<const std::int64_t> positions) const;
    // Sorted unique atoms (asunique).
    AtomGroup unique() const;
    // Sorted unique residue indices of the atoms (ag.residues).
    std::vector<std::int32_t> residues() const;
    // Bonds with at least one atom in the group (ag.bonds). Throws NoDataError without bonds.
    std::vector<Bond> bonds() const;

    bool operator==(const AtomGroup& other) const {
        return m_topology == other.m_topology && m_indices == other.m_indices;
    }

private:
    std::shared_ptr<const Topology> m_topology;
    std::vector<std::int32_t> m_indices;
};

// a | b: sorted unique union. Both must belong to the same topology (std::invalid_argument).
AtomGroup operator|(const AtomGroup& a, const AtomGroup& b);

// The atoms of residues, residue by residue in the given order (ResidueGroup.atoms).
AtomGroup residueAtoms(const std::shared_ptr<const Topology>& topology, std::span<const std::int32_t> residues);

// Per residue (first atom of each): the chain index MDATK colours and draws chains by
// (analysis/structure.py chain_ids) and the chain name (chain_names).
std::vector<std::int64_t> chainIds(const Topology& topology, std::span<const std::int32_t> residues);
std::vector<std::string> chainNames(const Topology& topology, std::span<const std::int32_t> residues);

// One frame: positions in Å (3 per atom, float32), the unit cell, the time in ps.
struct Frame {
    std::int64_t index = 0;  // 0-based frame of the (concatenated) trajectory
    double time = 0.0;  // ps: the file's time when stored (float32 widened), else index * dt
    std::vector<float> positions;
    std::optional<std::array<float, 6>> box;  // a, b, c (Å), alpha, beta, gamma (degrees)
};

// Frames of a system, with a current frame as MDAnalysis' trajectory has one. Implementations
// read files (io) or hold frames in memory.
class Trajectory {
public:
    virtual ~Trajectory() = default;
    virtual std::int64_t frameCount() const = 0;
    // ps between frames: the file's time step, or the one given when it was opened
    virtual double dt() const = 0;
    // Makes frame `index` (0-based; negative counts from the end) current and returns it.
    // Throws std::out_of_range.
    virtual const Frame& seek(std::int64_t index) = 0;
    virtual const Frame& current() const = 0;
};

// Frames held in memory (MDAnalysis' MemoryReader); times default to index * dt.
class MemoryTrajectory : public Trajectory {
public:
    // `frames`: positions (3 per atom) and boxes; index and, unless `keepTimes`, time are set here.
    MemoryTrajectory(std::vector<Frame> frames, double dt = 1.0, bool keepTimes = false);
    std::int64_t frameCount() const override;
    double dt() const override { return m_dt; }
    const Frame& seek(std::int64_t index) override;
    const Frame& current() const override;

private:
    std::vector<Frame> m_frames;
    double m_dt;
    std::size_t m_current = 0;
};

// A loaded system (MDAnalysis Universe): topology and, when coordinates were loaded, frames.
struct System {
    std::shared_ptr<const Topology> topology;
    std::unique_ptr<Trajectory> trajectory;

    AtomGroup atoms() const { return AtomGroup(topology); }
    // The current frame, or nullptr without coordinates.
    const Frame* currentFrame() const { return trajectory ? &trajectory->current() : nullptr; }
};

}  // namespace mdbridge::core
