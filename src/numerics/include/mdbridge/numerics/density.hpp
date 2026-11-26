#pragma once
// Density-map helpers of mdatk/density.py (the owner's code) that the window and the density fit share:
// density_centroid. Moved here from render/panels (milestone 6) so that chem may use it.

#include "mdbridge/core/rigid_transform.hpp"  // Vec3
#include "mdbridge/io/mrc.hpp"                // DensityMap

namespace mdbridge::numerics {

// density_centroid(map, level): the centre (Å) of the map density above `level` (voxel weights
// value - level, non-finite voxels carry none); the centre of the map box when nothing is above it.
// With the model centre it gives MDATK's centre_transform (the window's *Centre on model*).
core::Vec3 densityCentroid(const io::DensityMap& map, double level);

}  // namespace mdbridge::numerics
