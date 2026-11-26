// MDATK's rigid-body fit of a model into a map (mdatk/builder.py v0.0.4.1 _atom_element,
// model_atoms, fit_to_map; MDATK's own code, ported) on chem::fitModel, with CH-10 (non-finite
// resolutions refused) and CH-11 (in fitModel).

#include "mdbridge/chem/fit_to_map.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/mrc.hpp"

namespace mdbridge::chem {

namespace {

using io::BuildError;

// re.sub(r"[\d+\-']", "", text) for ASCII digits (atom and residue names are ASCII)
std::string withoutDigitsSigns(std::string_view text, bool primes) {
    std::string out;
    for (const char c : text) {
        if ((c >= '0' && c <= '9') || c == '+' || c == '-' || (primes && c == '\'')) {
            continue;
        }
        out.push_back(c);
    }
    return out;
}

std::string format(const char* spec, double x) {
    char text[64];
    std::snprintf(text, sizeof text, spec, x);
    return text;
}

// f"{x:g}"
std::string general(double x) {
    return format("%g", x);
}

}  // namespace

std::string atomElement(const io::PdbAtom& atom, std::string_view resname) {
    if (!atom.element.empty()) {
        return atom.element;
    }
    const std::string bare = io::pythonUpper(withoutDigitsSigns(atom.name, true));
    if (core::decodeUtf8(bare).size() == 2 && bare == io::pythonUpper(withoutDigitsSigns(resname, false))) {
        return bare;
    }
    return core::encodeUtf8(core::decodeUtf8(bare).substr(0, 1));
}

ModelAtoms modelAtoms(const io::PdbStructure& structure, bool heavyOnly) {
    ModelAtoms out;
    for (const io::PdbResidue* res : structure.residues()) {
        for (const io::PdbAtom& atom : res->atoms) {
            if (!(atom.alt.empty() || atom.alt == "A" || atom.alt == "1")) {
                continue;
            }
            const std::string el = io::pythonUpper(atomElement(atom, res->name));
            if (heavyOnly && (el == "H" || el == "D")) {
                continue;
            }
            out.xyz.push_back(atom.xyz());
            out.elements.push_back(el);
        }
    }
    return out;
}

FitToMapResult fitToMap(const PdbSource& source, const std::filesystem::path& mapPath, double resolution,
                        const FitToMapOptions& options) {
    if (!std::isfinite(resolution)) {  // CH-10
        throw BuildError("the resolution must be a positive number of Å");
    }
    if (std::find(kFitSearches.begin(), kFitSearches.end(), options.search) == kFitSearches.end()) {
        throw BuildError("search must be one of global, local");
    }
    if (!(resolution > 0)) {
        throw BuildError("enter the map resolution (Å)");
    }
    const std::string text = source.content();
    io::PdbStructure s = io::readPdbText(text);
    const ModelAtoms atoms = modelAtoms(s, options.heavyOnly);
    if (atoms.xyz.size() < 3) {
        throw BuildError("the model has fewer than three atoms to fit");
    }
    std::shared_ptr<const io::DensityMap> map;
    try {
        map = io::loadMapCached(mapPath);
    } catch (const std::filesystem::filesystem_error& e) {
        throw BuildError("density map: " + core::pythonOsErrorText(e));
    } catch (const io::MapError& e) {
        throw BuildError(std::string("density map: ") + e.what());
    }
    FitToMapResult out;
    out.notes = s.notes;
    core::Vec3 centre{};  // xyz.mean(axis=0)
    for (const core::Vec3& p : atoms.xyz) {
        for (std::size_t c = 0; c < 3; ++c) {
            centre[c] += p[c];
        }
    }
    for (double& c : centre) {
        c /= static_cast<double>(atoms.xyz.size());
    }
    const core::Vec3 extent = map->extent();
    bool inside = true;
    for (std::size_t c = 0; c < 3; ++c) {
        inside = inside && centre[c] >= map->origin[c] && centre[c] <= extent[c];
    }
    const bool global = options.search == "global";
    if (!global && !inside) {
        out.notes.emplace_back(
            "WARNING: the model centre lies outside the map box; a local fit needs a rough placement first (use "
            "the global search)");
    }
    try {
        out.fit = fitModel(*map, atoms.xyz, atoms.elements, resolution, std::nullopt,
                           global ? FitSearch::Global : FitSearch::Manual, options.angularStep, options.control);
    } catch (const io::FitError& e) {
        throw BuildError(std::string("density fit: ") + e.what());
    }
    const core::RigidTransform& transform = out.fit.transform;
    std::optional<std::string> moved = io::transformPdbText(text, transform.rotation, transform.translation);
    if (!moved) {  // columns damaged (e.g. ChimeraX): rewrite through the PDB tools
        for (std::vector<io::PdbResidue>& seg : s.segments) {
            for (io::PdbResidue& r : seg) {
                for (io::PdbAtom& a : r.atoms) {
                    const core::Vec3 x = transform.apply(a.xyz());
                    a.tail = a.tail.substr(0, 3) + io::newAtomTail(x).substr(3, 24) + a.tail.substr(27);
                    a.anisou.clear();
                }
            }
        }
        moved = chem::writePdb(s, false);
        out.notes.emplace_back(
            "Rewrote the file in standard PDB columns (the input had damaged ATOM records); ANISOU records dropped");
    }
    out.text = std::move(*moved);
    const core::Vec3 c2 = transform.apply(centre);
    const double dx = c2[0] - centre[0];
    const double dy = c2[1] - centre[1];
    const double dz = c2[2] - centre[2];
    double sq = dx * dx;
    sq += dy * dy;
    sq += dz * dz;
    const double shift = std::sqrt(sq);
    out.cc = out.fit.cc;
    const std::string label = global ? "global search" : "local refinement";
    out.notes.push_back("Density fit (" + label + ") into " + core::pathToUtf8(mapPath.filename()) + " at " +
                        general(resolution) + " Å, " + std::to_string(atoms.xyz.size()) + " atoms: CC " +
                        format("%.3f", out.fit.cc) + (out.fit.message.empty() ? "" : " (" + out.fit.message + ")"));
    out.notes.push_back("Model moved by a " + format("%.1f", core::rotationAngle(transform.rotation)) +
                        "° rotation and a " + format("%.1f", shift) + " Å shift of its centre");
    for (const std::string& line : core::pythonSplitlines(text)) {
        if (line.rfind("CRYST1", 0) == 0) {
            out.notes.emplace_back("CRYST1 / SCALE records kept as they were (they may not describe the map's cell)");
            break;
        }
    }
    return out;
}

}  // namespace mdbridge::chem
