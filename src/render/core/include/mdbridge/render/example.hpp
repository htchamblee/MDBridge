#pragma once
// Placeholder data for laying out the figure before anything is computed (MDATK's
// mdatk/render/example.py): the GUI's layout preview draws the real panels with this example
// dataset: smooth sine-like series (seeded with zlib.crc32 of the data key, numpy's PCG64 ported
// bit for bit, numpy_random.hpp), residue heatmaps and maps sized from the residue groups, and
// either the real reference-frame trace (when one has been loaded) or a schematic duplex. Nothing
// touches a trajectory. Qt-free.
//
// Differences from MDATK (docs/DIFFERENCES.md): a structure with a single strand gets no example
// base pairs instead of an IndexError (RN-21); residue groups are expanded only up to the 400
// residues used and the schematic sequence repeats for groups of more than 120 residues
// (MDATK: memory for every number of a huge range, IndexError beyond 120; RN-25).

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/numpy_random.hpp"
#include "mdbridge/render/orientation.hpp"

namespace mdbridge::render {

inline constexpr int kExampleFrames = 120;  // N_FRAMES
inline constexpr std::size_t kExampleGroupLimit = 400;  // residues per group used

// _rng(key): default_rng(zlib.crc32(key.encode("utf-8"))).
NumpyGenerator exampleRng(std::string_view key);

// group_resids(project, name, fallback): residue numbers of the first group named `name` with
// residues (ranges in order, at most 400; chains ignored); 1..fallback for a missing group, one
// without residues or with a spec that cannot be read.
std::vector<std::int64_t> groupResids(const core::Project& project, std::string_view name, int fallback = 12);

// The per-residue part of a structure: what Pipeline.structure_preview() returns for the reference
// frame (the GUI's "Load real frame"), or a schematic trace.
struct TraceSource {
    std::vector<Vec3d> coords;  // one frame, per residue
    std::vector<std::int64_t> chain;
    std::vector<std::int64_t> groupIndex;
    std::vector<std::int64_t> resid;
    std::vector<std::array<std::int64_t, 2>> links;
    std::vector<Vec3d> ions;
    std::optional<std::vector<std::string>> resname;
    std::optional<std::vector<std::string>> chainName;
    std::optional<std::vector<std::vector<bool>>> groupMember;  // per group, per residue
    std::vector<std::string> labels;                            // residue labels ("" none)
};

struct AtomsSource {
    std::vector<Vec3d> coords;  // one frame, per atom
    std::vector<std::string> names;
    std::vector<std::string> elements;
    std::vector<std::string> resnames;
    std::vector<std::int64_t> row;  // trace row of each atom
    std::vector<std::array<std::int64_t, 2>> bonds;
};

struct StructurePreview {
    TraceSource trace;
    std::optional<AtomsSource> atoms;
};

// schematic_trace: a two-strand helical trace of the first two shown residue groups
// (visual.structure_groups, else the first two groups, else "A").
TraceSource schematicTrace(const core::Project& project);
// schematic_structure: a two-strand duplex with heavy atoms (the demo's nucleotide geometry), at
// most 40 residues per strand: a C4' trace and the atoms.
StructurePreview schematicStructure(const core::Project& project);

// example_dataset(project, preview, n_frames): placeholder values for every data key the layout's
// panels use; with a structure, sequence or secondary panel also trace, basepairs, atoms (with
// atoms) and proximity results. nFrames >= 1 (core::ValueError otherwise).
core::Dataset exampleDataset(const core::Project& project, const StructurePreview* preview = nullptr,
                             int nFrames = kExampleFrames);

}  // namespace mdbridge::render
