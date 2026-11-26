#include "project_view.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <optional>
#include <utility>

#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/render/panels/structure.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;

namespace {

using core::AnalysisResult;
using core::NumArray;

const NumArray& required(const AnalysisResult& result, std::string_view name) {
    const NumArray* a = result.array(name);
    if (a == nullptr) {
        throw core::ValueError("the " + result.name + " result has no " + std::string(name) + " array");
    }
    return *a;
}

std::vector<std::int64_t> ints(const NumArray& a) {
    const std::vector<double> d = a.toDouble();
    std::vector<std::int64_t> out(d.size());
    std::transform(d.begin(), d.end(), out.begin(), [](double v) { return static_cast<std::int64_t>(v); });
    return out;
}

// the first frame of a frames x rows x 3 array
std::vector<render::Vec3d> firstFrame(const NumArray& a) {
    const std::vector<std::size_t>& shape = a.shape();
    const std::size_t rows = shape.size() == 3 && shape[0] > 0 ? shape[1] : 0;
    const std::vector<double> d = a.toDouble();
    std::vector<render::Vec3d> out(rows);
    for (std::size_t i = 0; i < rows; ++i) out[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
    return out;
}

std::vector<std::array<std::int64_t, 2>> pairs(const NumArray& a) {
    const std::vector<std::int64_t> v = ints(a);
    std::vector<std::array<std::int64_t, 2>> out;
    for (std::size_t i = 0; i + 1 < v.size(); i += 2) out.push_back({v[i], v[i + 1]});
    return out;
}

std::vector<std::string> strings(const NumArray& a) {
    if (a.dtype() == NumArray::DType::String) {
        const std::span<const std::string> s = a.values<std::string>();
        return {s.begin(), s.end()};
    }
    std::vector<std::string> out;
    for (const std::int64_t v : ints(a)) out.push_back(std::to_string(v));
    return out;
}

std::vector<std::string> metaStrings(const core::Table& meta, std::string_view key) {
    std::vector<std::string> out;
    for (const auto& [k, v] : meta) {
        if (k == key) {
            if (const core::Array* a = v.get<core::Array>()) {
                for (const core::Value& e : *a) out.push_back(e.stringOr());
            }
        }
    }
    return out;
}

}  // namespace

core::Project structureViewProject(const core::Project& project) {
    core::Project p = project;
    p.layout.preset = "single";
    p.layout.panels = {core::PanelSpec{"structure", {}, "", {}}};
    if (p.visual.representations.empty()) {
        core::RepSpec rep;
        rep.style = std::string(kDefaultViewStyle);
        p.visual.representations.push_back(rep);
    }
    p.visual.renderer = "software";
    return p;
}

std::vector<std::string> referenceAnalyses(const core::Project& project) {
    return analysis::panelRequiredAnalyses(structureViewProject(project));
}

ReferenceFrame loadReferenceFrame(const core::Project& project, const core::Reporter& reporter) {
    ReferenceFrame out;
    const core::InputParams& in = project.inputs;
    std::string topology = in.topology;
    std::vector<std::string> trajectories = in.trajectories;
    {
        // structure_preview: cleaned files when a matching result exists, else the raw files (no cpptraj run)
        const analysis::Pipeline pipeline(project, core::Reporter());
        const analysis::ResolvedInputs resolved = pipeline.cleanedIfAvailable();
        topology = resolved.topology;
        trajectories = resolved.trajectories;
        out.cleaned = !resolved.cleaningKey.empty();
        if (project.cleaning.enabled && !out.cleaned) {
            out.notes.push_back("no cleaned trajectory of these inputs yet: the structure is shown from the raw files");
        }
    }
    reporter.checkCancel();
    std::vector<fs::path> trajectoryPaths;
    for (const std::string& t : trajectories) trajectoryPaths.push_back(core::pathFromUtf8(t));
    core::System system = io::loadSystem(core::pathFromUtf8(topology), trajectoryPaths, in.dt_ps);
    for (const std::string& w : system.warnings) {
        out.notes.push_back(w);
        reporter.log("WARNING: " + w);
    }
    if (!system.trajectory) {
        throw core::ValueError("No trajectory loaded into Universe");
    }
    const std::vector<std::int64_t> idx =
        io::frameIndices(system.trajectory->frameCount(), in.start, in.stop, in.step);
    const std::vector<std::int64_t> first{idx.front()};
    out.frame = idx.front();
    out.topology = topology;

    const std::vector<std::string> names = referenceAnalyses(project);
    std::vector<std::unique_ptr<analysis::AnalysisContext>> contexts;
    std::vector<std::unique_ptr<analysis::Analysis>> analyses;
    for (const std::string& n : names) {
        contexts.push_back(std::make_unique<analysis::AnalysisContext>(system, analysis::instanceProject(project, n),
                                                                        topology, reporter));
        analyses.push_back(analysis::analysisFor(n).create(*contexts.back(), n));
    }
    analysis::RunResult run;
    if (!analyses.empty()) {
        run = analysis::runAnalyses(*contexts.front(), analyses, first, &reporter, analysis::RunOptions{1});
    } else {
        run.timePs.push_back(system.trajectory->seek(out.frame).time);
    }
    out.dataset = core::Dataset(first, run.timePs, std::move(run.results), in.time_unit);
    return out;
}

render::StructurePreview structurePreview(const core::Dataset& reference) {
    const AnalysisResult* tr = reference.result("trace");
    if (tr == nullptr) {
        throw core::ValueError("the reference frame has no trace");
    }
    render::StructurePreview p;
    render::TraceSource& t = p.trace;
    t.coords = firstFrame(required(*tr, "coords"));
    if (const NumArray* ions = tr->array("ion_coords")) t.ions = firstFrame(*ions);
    t.links = pairs(required(*tr, "links"));
    t.chain = ints(required(*tr, "chain"));
    t.groupIndex = ints(required(*tr, "group_index"));
    t.resid = ints(required(*tr, "resid"));
    if (const NumArray* a = tr->array("resname")) t.resname = strings(*a);
    if (const NumArray* a = tr->array("chain_name")) t.chainName = strings(*a);
    if (const NumArray* gm = tr->array("group_member")) {
        const std::vector<double> m = gm->toDouble();
        const std::size_t groups = gm->shape().size() == 2 ? gm->shape()[0] : 0;
        const std::size_t res = t.coords.size();
        std::vector<std::vector<bool>> member(groups, std::vector<bool>(res));
        for (std::size_t g = 0; g < groups; ++g) {
            for (std::size_t r = 0; r < res && g * res + r < m.size(); ++r) member[g][r] = m[g * res + r] != 0.0;
        }
        t.groupMember = std::move(member);
    }
    t.labels = metaStrings(tr->meta, "residue_labels");
    if (const AnalysisResult* at = reference.result("atoms")) {
        render::AtomsSource& s = p.atoms.emplace();
        s.coords = firstFrame(required(*at, "coords"));
        s.names = metaStrings(at->meta, "atom_names");
        s.elements = metaStrings(at->meta, "elements");
        s.resnames = metaStrings(at->meta, "resnames");
        s.row = ints(required(*at, "row"));
        s.bonds = pairs(required(*at, "bonds"));
    }
    return p;
}

StructureView structureView(const core::Project& project, const core::Dataset& reference) {
    const core::Project viewProject = structureViewProject(project);
    render::panels::Context ctx(reference, viewProject);
    render::Figure fig(400, 400, 72.0);
    const core::PanelSpec& spec = viewProject.layout.panels.front();
    const std::unique_ptr<render::panels::Panel> panel =
        render::panels::Panel::make<render::panels::StructurePanel>(fig, render::SubplotSpec{}, spec, ctx);
    auto& structure = static_cast<render::panels::StructurePanel&>(*panel);
    if (!structure.ok()) {
        throw core::ValueError("the structure view has no data: the reference frame lacks the trace or the atoms");
    }
    std::optional<render::panels::StructurePanel::StageFrame> stage = structure.stageFrame(0);
    if (!stage || stage->items.empty()) {
        throw core::ValueError("no 3D representation could be drawn for this structure");
    }
    StructureView view;
    view.items = std::move(stage->items);
    view.viewport = stage->viewport;
    view.shading = std::move(stage->shading);
    for (std::size_t k = 0; k < 3; ++k) view.background[k] = static_cast<float>(view.shading.fog_rgb[k]);
    if (const AnalysisResult* tr = reference.result("trace")) {
        if (const NumArray* c = tr->array("coords"); c && c->shape().size() == 3) view.residues = c->shape()[1];
    }
    if (const AnalysisResult* at = reference.result("atoms")) {
        if (const NumArray* c = at->array("coords"); c && c->shape().size() == 3) view.atoms = c->shape()[1];
    }
    return view;
}

std::string previewWatermark(bool referenceFrame) {
    return referenceFrame ? "Example data – layout preview (structure: reference frame)"
                          : "Example data – layout preview";
}

PreviewSize previewSize(const core::Project& project, int width, int height) {
    const core::VideoParams& vp = project.output.video;
    const double vw = static_cast<double>(std::max<std::int64_t>(vp.width, 1));
    const double vh = static_cast<double>(std::max<std::int64_t>(vp.height, 1));
    const double scale = std::min(width / vw, height / vh);
    PreviewSize size;
    size.width = std::max(2, static_cast<int>(vw * scale));
    size.height = std::max(2, static_cast<int>(vh * scale));
    size.dpi = static_cast<double>(vp.dpi) * scale;
    return size;
}

std::size_t previewPosition(std::size_t frames) {
    return static_cast<std::size_t>(kPreviewFrame * static_cast<double>(frames > 0 ? frames - 1 : 0));
}

render::panels::BuiltFigure previewFigure(const core::Project& project, const core::Dataset& example,
                                          const render::FontSet& fonts, int width, int height,
                                          const std::string& watermark) {
    const PreviewSize size = previewSize(project, width, height);
    render::panels::BuiltFigure built =
        render::panels::buildFigure(example, project, fonts, size.width, size.height, size.dpi, watermark);
    if (example.frameCount() > 0) {
        built.update(previewPosition(example.frameCount()), fonts);
    }
    return built;
}

}  // namespace mdbridge::app
