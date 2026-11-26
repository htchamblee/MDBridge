#pragma once
// MDATK's structure builder (mdatk/builder.py v0.0.4.1, MDATK's own code, ported), part 4: the
// rigid-body fit of a model into a density map (`build fit`): the model's heavy atoms (first
// alternate location) fitted with chem::fitModel (global search or local refinement), then only the
// coordinates of the input text changed (io::transformPdbText: ATOM / HETATM coordinates, ANISOU
// rotated, every other line verbatim). An input whose atom records lack fixed-column coordinates
// (ChimeraX-damaged lines the PDB reader repairs) is rewritten in standard columns instead.
// Differences from MDATK v0.0.4.1 (docs/DIFFERENCES.md):
//   CH-10 a non-finite resolution: BuildError "the resolution must be a positive number of Å"
//         (MDATK failed in NumPy)
//   CH-11 non-finite map voxels carry no density (chem::fitModel), the fit message ends "; N
//         non-finite map voxel(s) read as no density"
// The fit itself is within D10's tolerance of MDATK's (chem/density_fit.hpp).
//
// Callers: `build fit` (milestone 6 W2-C), the builder window's *Fit to map* page (W3-B).

#include <array>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/chem/density_fit.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

// builder.SEARCHES
inline constexpr std::array<std::string_view, 2> kFitSearches = {"global", "local"};

// _atom_element: the element column, else the name without digits, + - and ' (an ion named after
// its residue keeps both letters, MG in MG; otherwise the first letter). Not upper-cased.
std::string atomElement(const io::PdbAtom& atom, std::string_view resname);

struct ModelAtoms {
    std::vector<core::Vec3> xyz;
    std::vector<std::string> elements;  // upper case
};
// model_atoms(structure, heavy_only): the atoms of the first alternate location ("", "A" or "1"),
// without H / D atoms when heavyOnly, in file order.
ModelAtoms modelAtoms(const io::PdbStructure& structure, bool heavyOnly = true);

struct FitToMapOptions {
    std::string search = "global";  // "global" (exhaustive search, then refinement) or "local"
    double angularStep = kDefaultAngularStep;  // degrees, the global search's rotation spacing
    bool heavyOnly = true;
    FitControl control;  // progress, cancel (core::Cancelled), threads
};

struct FitToMapResult {
    std::string text;  // the input text with the moved coordinates (or the rewritten file)
    double cc = 0.0;
    Notes notes;       // read notes, the outside-the-box WARNING, "Density fit (...) ...", "Model moved by ..."
    FitResult fit;     // fitModel's result (transform: world -> map)
};

// fit_to_map(source, map_path, resolution, search, angular_step, heavy_only). BuildError: CH-10's
// message, "search must be one of global, local", "enter the map resolution (Å)" (zero, negative),
// "the model has fewer than three atoms to fit", "density map: [Errno 2] No such file or directory:
// 'x.mrc'" (and the MapError messages), "density fit: ..." (FitError: "the angular step must be at
// least 3 degrees (got 2.0)" ...), "the fitted coordinates do not fit the PDB columns"; the
// io::PdbError of reading the model; core::Cancelled.
FitToMapResult fitToMap(const PdbSource& source, const std::filesystem::path& mapPath, double resolution,
                        const FitToMapOptions& options = {});

}  // namespace mdbridge::chem
