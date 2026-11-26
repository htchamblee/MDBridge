#pragma once
// MDATK's synthetic demonstration system (mdatk/demo.py, MDATK's own code, ported): an idealised
// 12-bp RNA "aptamer / target" duplex with three Mg2+ ions, Na+ counter-ions and a water box,
// written as an Amber topology (demo.parm7), a NetCDF trajectory (demo.nc) and the project file
// demo.mdatk.toml. The trajectory has thermal noise, a slow rigid-body drift, a transient fraying
// of one duplex end and Mg ions that bind, bridge and leave phosphates. The files are MDATK
// v0.0.4.1's byte for byte for the same arguments (numpy's default_rng(seed) through
// render::NumpyGenerator; the arithmetic in numpy's order, without fused multiply-adds), apart
// from the project file's `saved` time.
//
// Used by the `demo` command and by the interface (Create demo system…, Demo project…). No Qt.

#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::app {

struct DemoOptions {
    std::int64_t frames = 120;  // --n-frames
    double dtPs = 10.0;         // --dt
    bool solvent = true;        // not --no-solvent
    std::string sequence = "GCAUGGCAUCGA";
    std::uint64_t seed = 7;
};

// build_demo's result: the paths as os.path.join(out_dir, name) gives them (relative when
// out_dir is), in MDATK's order: topology, trajectory, project, preset (the project again).
struct DemoFiles {
    std::string topology;
    std::string trajectory;
    std::string project;
    std::vector<std::pair<std::string, std::string>> entries() const {
        return {{"topology", topology}, {"trajectory", trajectory}, {"project", project}, {"preset", project}};
    }
};

// demo_project(topology, trajectory, out_dir, n_bp, cleaning, mid_time_ps): the demo's settings
// (groups 1-n / n+1-2n, screened energy, side view, PDB frames "0, last" and the middle time,
// layout large_plus_4 with an Mg occupancy heatmap).
core::Project demoProject(const std::string& topology, const std::string& trajectory, const std::string& outDir,
                          std::int64_t basePairs, bool cleaning, double midTimePs = 300.0);

// build_demo(out_dir, ...): creates out_dir and writes the three files (the project with
// cleaning on when there is solvent and cpptraj is on PATH). `outDir` is UTF-8 text, joined as
// os.path.join does for the returned paths. `reporter`: progress of the frames and its cancel
// flag (core::Cancelled; nothing is written before the frames are done). Throws core::ValueError,
// before anything is written, for fewer than one frame, a time step that is not a positive finite
// number (MDATK: a crash after writing the topology, or files with times of 0 / negative / nan
// ps; DIFFERENCES CMD-01) and a sequence of fewer than four of A, C, G, U;
// std::filesystem::filesystem_error.
DemoFiles buildDemo(const std::string& outDir, const DemoOptions& options = {},
                    const core::Reporter& reporter = core::Reporter());

}  // namespace mdbridge::app
