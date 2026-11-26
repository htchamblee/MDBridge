#include "density_pipeline.hpp"

#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/molecule/density_shadow.hpp"
#include "mdbridge/render/panels/structure.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;
using analysis::PipelineError;

namespace {

// f"{x:g}"
std::string formatG(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%g", x);
    return text;
}

std::string fixed3(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%.3f", x);
    return text;
}

// frames x rows x 3 float array: the rows of frame `k` in double
std::vector<core::Vec3> frameRows(const core::NumArray& a, std::size_t k) {
    const std::vector<std::size_t>& shape = a.shape();
    if (shape.size() != 3 || shape[0] <= k || shape[2] != 3) {
        throw std::logic_error("coordinates of an unexpected shape");
    }
    const std::vector<double> d = a.toDouble();
    std::vector<core::Vec3> out(shape[1]);
    const std::size_t base = k * shape[1] * 3;
    for (std::size_t i = 0; i < shape[1]; ++i) {
        out[i] = {d[base + 3 * i], d[base + 3 * i + 1], d[base + 3 * i + 2]};
    }
    return out;
}

std::string searchLabel(const std::string& search) { return search == "global" ? "global search" : "local refinement"; }

}  // namespace

DensityFrameAtoms densityFrameAtoms(analysis::Pipeline& pipeline, std::string_view frame,
                                    const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : pipeline.reporter();
    const core::Project& p = pipeline.project();
    core::System& system = pipeline.system();
    const std::int64_t frameCount = system.trajectory->frameCount();
    const std::vector<std::int64_t> idx = io::frameIndices(frameCount, p.inputs.start, p.inputs.stop, p.inputs.step);
    if (idx.empty()) {
        throw PipelineError("no frames selected");
    }
    const std::int64_t ref = idx.front();
    const std::string spec = core::pythonLower(core::pythonStrip(frame.empty() ? std::string_view("ref") : frame));
    const std::int64_t f = spec.empty() || spec == "ref" ? ref : core::pythonInt(spec);
    if (f < 0 || f >= frameCount) {
        throw PipelineError("density fit frame " + std::to_string(f) + " is outside the trajectory (0-" +
                            std::to_string(frameCount - 1) + ")");
    }
    const std::vector<std::int64_t> frames = f == ref ? std::vector<std::int64_t>{ref} : std::vector<std::int64_t>{ref, f};
    const analysis::ResolvedInputs& inputs = pipeline.resolveInputs();
    analysis::AnalysisContext context(system, p, inputs.topology, rep);
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    analyses.push_back(analysis::analysisFor("atoms").create(context, "atoms"));
    analyses.push_back(analysis::analysisFor("trace").create(context, "trace"));
    const core::Reporter silent;
    const analysis::RunResult run = analysis::runAnalyses(context, analyses, frames, &silent);
    const core::AnalysisResult& atoms = run.results.at(0);
    const core::AnalysisResult& trace = run.results.at(1);
    DensityFrameAtoms out;
    const core::NumArray* coords = atoms.array("coords");
    const core::NumArray* traceCoords = trace.array("coords");
    if (coords == nullptr || traceCoords == nullptr) {
        throw std::logic_error("the atoms and trace analyses gave no coordinates");
    }
    out.xyz = frameRows(*coords, frames.size() - 1);
    for (const auto& [key, value] : atoms.meta) {
        if (key == "elements") {
            if (const core::Array* list = value.get<core::Array>()) {
                for (const core::Value& e : *list) out.elements.push_back(e.stringOr());
            }
        }
    }
    out.frame = f;
    out.reference = frameRows(*traceCoords, 0);
    return out;
}

chem::FitResult fitDensity(analysis::Pipeline& pipeline, std::optional<std::string> search,
                           const std::optional<std::vector<double>>& start, const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : pipeline.reporter();
    core::DensityParams d = pipeline.project().density;
    if (d.map_path.empty()) {
        throw PipelineError("no density map selected");
    }
    std::shared_ptr<const io::DensityMap> map;
    try {
        map = io::loadMapCached(core::pathFromUtf8(d.map_path));
    } catch (const fs::filesystem_error& error) {
        throw PipelineError("density map: " + core::pythonOsErrorText(error));
    } catch (const io::MapError& error) {
        throw PipelineError(std::string("density map: ") + error.what());
    }
    const std::string how = search && !search->empty() ? *search : d.search;
    rep.log("Density fit (" + searchLabel(how) + ") of frame " + (d.fit_frame.empty() ? "ref" : d.fit_frame) +
            " into " + core::pathBasename(d.map_path) + " at " + formatG(d.resolution) + " Å");
    const DensityFrameAtoms atoms = densityFrameAtoms(pipeline, d.fit_frame, &rep);
    std::optional<core::RigidTransform> t0;
    if (how != "global") {
        if (start) {
            t0 = core::RigidTransform::fromList(*start);
        } else if (d.transform.size() == 12) {
            try {
                const render::molecule::AnchoredPlacement placed = render::molecule::anchored(d.transform, d.anchor, atoms.reference);
                if (!placed.note.empty()) {
                    rep.log("WARNING: " + placed.note);
                }
                t0 = placed.transform;
            } catch (const core::ValueError& error) {
                throw PipelineError(std::string("density map: stored placement: ") + error.what());
            }
        } else {
            core::Vec3 centre{0, 0, 0};
            for (const core::Vec3& x : atoms.xyz) {
                for (int k = 0; k < 3; ++k) centre[static_cast<std::size_t>(k)] += x[static_cast<std::size_t>(k)];
            }
            for (double& c : centre) c /= static_cast<double>(atoms.xyz.empty() ? 1 : atoms.xyz.size());
            t0 = chem::centreTransform(*map, render::panels::finiteMapLevel(*map, d.level), centre);  // CH-11
        }
    }
    chem::FitResult result;
    try {
        chem::FitControl control;
        control.reporter = &rep;
        result = chem::fitModel(*map, atoms.xyz, atoms.elements, d.resolution, t0,
                                how == "global" ? chem::FitSearch::Global : chem::FitSearch::Manual, d.angular_step,
                                control);
    } catch (const chem::FitError& error) {
        throw PipelineError(std::string("density fit: ") + error.what());
    }
    d.transform = result.transform.toList();
    d.anchor.clear();
    for (const core::Vec3& a : core::anchorPoints(atoms.reference)) {
        d.anchor.push_back({a[0], a[1], a[2]});
    }
    d.fit_cc = core::pythonRound(result.cc, 4);
    d.fit_info = searchLabel(how) + ", frame " + std::to_string(atoms.frame) + ", " + std::to_string(atoms.xyz.size()) +
                 " atoms, " + formatG(d.resolution) + " Å: CC " + fixed3(result.cc) +
                 (result.message.empty() ? std::string() : " (" + result.message + ")");
    rep.log("Density fit: " + d.fit_info);
    pipeline.setDensity(std::move(d));
    return result;
}

void ensureDensityFit(analysis::Pipeline& pipeline, const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : pipeline.reporter();
    const core::Project& p = pipeline.project();
    if (!core::densityActive(p) || p.density.transform.size() == 12) {
        return;
    }
    if (!p.density.auto_fit) {
        rep.log("WARNING: the density map has no placement and automatic fitting is off; it is not drawn");
        return;
    }
    try {
        fitDensity(pipeline, std::nullopt, std::nullopt, reporter);
    } catch (const core::Cancelled&) {
        throw;
    } catch (const PipelineError& error) {
        rep.log(std::string("WARNING: density map not drawn: ") + error.what());
    } catch (const fs::filesystem_error& error) {
        rep.log("WARNING: density map not drawn: " + core::pythonOsErrorText(error));
    } catch (const core::ValueError& error) {
        rep.log(std::string("WARNING: density map not drawn: ") + error.what());
    }
}

}  // namespace mdbridge::app
