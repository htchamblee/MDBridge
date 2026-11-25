// The structure view (MDATK's mdatk/render/panels.py StructurePanel, molecule/stage.py
// MoleculeStage / PixelImage / Viewport, molecule/renderers.py get_renderer; MDATK's own code,
// ported). Improvements over MDATK: structure.hpp.

#include "mdbridge/render/panels/structure.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/render/molecule/density_shadow.hpp"

namespace mdbridge::render::panels {

using core::AnalysisResult;
using core::NumArray;

StageRenderer::~StageRenderer() = default;
StageRendererFactory::~StageRendererFactory() = default;

namespace {
std::mutex g_factoryMutex;
std::shared_ptr<const StageRendererFactory> g_factory;
}  // namespace

void setStageRendererFactory(std::shared_ptr<const StageRendererFactory> factory) {
    const std::lock_guard<std::mutex> lock(g_factoryMutex);
    g_factory = std::move(factory);
}

std::shared_ptr<const StageRendererFactory> stageRendererFactory() {
    const std::lock_guard<std::mutex> lock(g_factoryMutex);
    return g_factory;
}

StageRendererChoice defaultStageRenderer(std::string_view name) {
    static const std::shared_ptr<const StageRenderer> software = std::make_shared<SoftwareStageRenderer>();
    if (name == "vtk" || name == "vulkan") {
        return {software, std::string(kNoGpuRendererNote)};  // RN-16
    }
    return {software, {}};
}

double finiteMapLevel(const io::DensityMap& map, double sigmas) {
    const bool allFinite = std::all_of(map.data.begin(), map.data.end(), [](float v) { return std::isfinite(v); });
    if (allFinite) {
        return io::mapLevel(map, sigmas);
    }
    // RN-33: mean and (population) standard deviation of the finite voxels, in double
    double sum = 0.0;
    std::size_t n = 0;
    for (const float v : map.data) {
        if (std::isfinite(v)) {
            sum += v;
            ++n;
        }
    }
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double mean = sum / static_cast<double>(n);
    double ss = 0.0;
    for (const float v : map.data) {
        if (std::isfinite(v)) {
            const double d = static_cast<double>(v) - mean;
            ss += d * d;
        }
    }
    return mean + sigmas * std::sqrt(ss / static_cast<double>(n));
}

std::string ionLegendLabel(const std::vector<std::string>& ionLabels) {
    static const std::map<std::string, std::string, std::less<>> names{
        {"MG", "Mg\xc2\xb2\xe2\x81\xba"}, {"NA", "Na\xe2\x81\xba"},       {"SOD", "Na\xe2\x81\xba"},
        {"K", "K\xe2\x81\xba"},           {"POT", "K\xe2\x81\xba"},        {"CA", "Ca\xc2\xb2\xe2\x81\xba"},
        {"CAL", "Ca\xc2\xb2\xe2\x81\xba"}, {"ZN", "Zn\xc2\xb2\xe2\x81\xba"}, {"MN", "Mn\xc2\xb2\xe2\x81\xba"},
        {"CL", "Cl\xe2\x81\xbb"},         {"CLA", "Cl\xe2\x81\xbb"},       {"LI", "Li\xe2\x81\xba"},
        {"RB", "Rb\xe2\x81\xba"},         {"CS", "Cs\xe2\x81\xba"},
    };
    std::set<std::string> found;
    for (const std::string& label : ionLabels) {
        std::string name = label;
        while (!name.empty() && (std::isdigit(static_cast<unsigned char>(name.back())) != 0)) name.pop_back();
        std::string key;
        for (const char c : name) {
            if (c != '+' && c != '-' && c != ' ') key.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
        }
        while (!key.empty() && (std::isdigit(static_cast<unsigned char>(key.back())) != 0)) key.pop_back();  // MG2
        const auto it = names.find(key);
        found.insert(it == names.end() ? std::string("Ions") : it->second);
    }
    return found.size() == 1 ? *found.begin() : std::string("Ions");
}

namespace {

const NumArray& required(const AnalysisResult& r, std::string_view name) {
    const NumArray* a = r.array(name);
    if (a == nullptr) {
        throw NoData(std::string(name));
    }
    return *a;
}

std::vector<std::int64_t> ints(const NumArray& a) {
    const std::vector<double> d = a.toDouble();
    std::vector<std::int64_t> out(d.size());
    std::transform(d.begin(), d.end(), out.begin(), [](double v) { return static_cast<std::int64_t>(v); });
    return out;
}

std::vector<Vec3d> vec3s(const std::vector<double>& d) {
    std::vector<Vec3d> out(d.size() / 3);
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
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

bool metaFlag(const core::Table& meta, std::string_view key) {
    for (const auto& [k, v] : meta) {
        if (k == key) {
            if (const bool* b = v.get<bool>()) return *b;
            if (const std::int64_t* i = v.get<std::int64_t>()) return *i != 0;
            return !v.isNull();
        }
    }
    return false;
}

bool finite(const Vec3d& p) { return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]); }

// np.linalg.norm of a 3-vector: sqrt of the squares summed in order
double norm3(const Vec3d& p) { return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); }

// x @ m.T for one point (m's rows dotted with x, k = 0, 1, 2 in order)
Vec3d rotate(const Vec3d& x, const Mat3& m) {
    Vec3d out{};
    for (std::size_t i = 0; i < 3; ++i) out[i] = x[0] * m[i][0] + x[1] * m[i][1] + x[2] * m[i][2];
    return out;
}

// np.argsort(keys, kind="stable")
std::vector<std::size_t> stableOrder(const std::vector<std::size_t>& rows, const std::vector<double>& keys) {
    std::vector<std::size_t> out = rows;
    std::stable_sort(out.begin(), out.end(), [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
    return out;
}

Rgba withAlphaOf(const std::array<double, 3>& rgb, double alpha) { return {rgb[0], rgb[1], rgb[2], alpha}; }

std::string trimmed(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n\f\v");
    if (b == std::string::npos) return {};
    return s.substr(b, s.find_last_not_of(" \t\r\n\f\v") - b + 1);
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string formatG(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%g", v);
    return buf;
}

// AxesImage of float RGBA data: the bytes to_rgba(bytes=True) shows ((x * 255) as uint8, float32)
std::vector<std::uint8_t> rgbaBytes(const std::vector<float>& rgba) {
    std::vector<std::uint8_t> out(rgba.size());
    for (std::size_t i = 0; i < rgba.size(); ++i) {
        const float x = rgba[i] * 255.0f;
        out[i] = static_cast<std::uint8_t>(x > 0.0f ? x : 0.0f);
    }
    return out;
}

// density_centroid(map, level): the centre of the density above the level (non-finite voxels carry
// none); the map box' centre when nothing is above it.
Vec3d densityCentroid(const io::DensityMap& map, double level) {
    const auto nx = static_cast<std::size_t>(map.shapeXyz[0]), ny = static_cast<std::size_t>(map.shapeXyz[1]),
               nz = static_cast<std::size_t>(map.shapeXyz[2]);
    double total = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                const double v = static_cast<double>(map.data[(z * ny + y) * nx + x]) - level;
                if (!(v > 0.0) || !std::isfinite(v)) continue;
                total += v;
                cx += v * static_cast<double>(x);
                cy += v * static_cast<double>(y);
                cz += v * static_cast<double>(z);
            }
        }
    }
    if (!(total > 0.0)) {
        return {map.origin[0] + 0.5 * static_cast<double>(map.shapeXyz[0] - 1) * map.voxel[0],
                map.origin[1] + 0.5 * static_cast<double>(map.shapeXyz[1] - 1) * map.voxel[1],
                map.origin[2] + 0.5 * static_cast<double>(map.shapeXyz[2] - 1) * map.voxel[2]};
    }
    return map.indexToXyz({cx / total, cy / total, cz / total});
}

}  // namespace

// ---- the 3D layer ---------------------------------------------------------------------------------

struct StructurePanel::Stage {
    std::unique_ptr<molecule::MoleculeScene> scene;
    std::vector<molecule::Vec3> atoms;  // frames x atoms, centred on the panel centre (float32 values)
    std::size_t nAtoms = 0;
    std::shared_ptr<const StageRenderer> renderer;
    std::optional<std::size_t> pos;
    std::vector<molecule::Colored> items;
    bool failed = false;
    std::size_t image = 0;  // artist index of the PixelImage
};

StructurePanel::StructurePanel(const core::PanelSpec& spec, Context& context)
    : Panel("structure", spec, context), m_renderers(stageRendererFactory()) {}

StructurePanel::~StructurePanel() = default;

std::string StructurePanel::stageRenderer() const { return m_stage ? m_stage->renderer->name() : std::string(); }
bool StructurePanel::stageFailed() const { return m_stage && m_stage->failed; }

std::span<const Vec3d> StructurePanel::view(std::size_t pos) const {
    return std::span<const Vec3d>(m_viewXyz).subspan(pos * m_res, m_res);
}
std::span<const Vec3d> StructurePanel::ionView(std::size_t pos) const {
    return std::span<const Vec3d>(m_ionViewXyz).subspan(pos * m_nIons, m_nIons);
}

void StructurePanel::setup(Figure& fig) {
    const core::Dataset& ds = ctx.dataset;
    const core::VisualParams& v = ctx.project.visual;
    const AnalysisResult* tr = ds.result("trace");
    if (tr == nullptr) {
        throw NoData("trace");
    }
    const NumArray& coordArray = required(*tr, "coords");
    const NumArray& ionArray = required(*tr, "ion_coords");
    const std::vector<std::size_t>& shape = coordArray.shape();
    m_frames = shape.empty() ? 0 : shape[0];
    m_res = shape.size() > 1 ? shape[1] : 0;
    m_nIons = ionArray.shape().size() > 1 ? ionArray.shape()[1] : 0;
    if (m_frames == 0) {
        throw NoData("no frames");  // RN-02: MDATK's IndexError
    }
    const std::vector<Vec3d> coords = vec3s(coordArray.toDouble());
    const std::vector<Vec3d> ions = vec3s(ionArray.toDouble());
    viewCentre(coords, m_frames, m_res);
    m_centred.resize(coords.size());
    m_finite.assign(coords.size(), true);
    for (std::size_t f = 0; f < m_frames; ++f) {
        for (std::size_t r = 0; r < m_res; ++r) {
            const std::size_t i = f * m_res + r;
            for (std::size_t c = 0; c < 3; ++c) m_centred[i][c] = coords[i][c] - m_center[c] - m_shift[f][c];
            m_finite[i] = finite(m_centred[i]);
        }
    }
    m_ionsCentred.resize(ions.size());
    std::vector<bool> ionFinite(ions.size(), true);
    for (std::size_t f = 0; f < m_frames && m_nIons > 0; ++f) {
        for (std::size_t r = 0; r < m_nIons; ++r) {
            const std::size_t i = f * m_nIons + r;
            for (std::size_t c = 0; c < 3; ++c) m_ionsCentred[i][c] = ions[i][c] - m_center[c] - m_shift[f][c];
            ionFinite[i] = finite(m_ionsCentred[i]);
        }
    }
    {
        const std::vector<std::int64_t> flat = ints(required(*tr, "links"));
        for (std::size_t i = 0; i + 1 < flat.size(); i += 2) m_links.push_back({flat[i], flat[i + 1]});
    }
    // Frame the model by its bounding sphere so any orientation fits without zooming; non-finite
    // points are left out (RN-03)
    double radius = -std::numeric_limits<double>::infinity();
    std::set<std::size_t> bad;
    for (std::size_t i = 0; i < m_centred.size(); ++i) {
        if (m_finite[i]) radius = std::max(radius, norm3(m_centred[i]));
        else bad.insert(i / std::max<std::size_t>(m_res, 1));
    }
    for (std::size_t i = 0; i < m_ionsCentred.size(); ++i) {
        if (ionFinite[i]) radius = std::max(radius, norm3(m_ionsCentred[i]));
        else bad.insert(i / std::max<std::size_t>(m_nIons, 1));
    }
    if (!bad.empty()) {
        std::string list;
        std::size_t k = 0;
        for (const std::size_t f : bad) {
            if (k == 10) break;
            list += (k++ ? ", " : "") + std::to_string(ds.frames().at(f));
        }
        if (bad.size() > 10) list += " and " + std::to_string(bad.size() - 10) + " more";
        warnOnce("structure view: non-finite coordinates in frame(s) " + list +
                 "; left out of the view (those residues are not drawn there)");
    }
    radius = std::max(radius, 1.0);
    m_radius = radius;
    const double zoom = v.zoom;
    const double half = radius / std::max(zoom, 0.05);  // MDATK's "view"; zoom > 1 shows a smaller window
    const double pad = 0.06 * half + 2.0 / std::max(zoom, 1.0);
    m_dataLo = {-half - pad, -half - pad};
    m_dataHi = {half + pad, half + pad};
    Axes& a = ax(fig);
    a.setXLim(m_dataLo[0], m_dataHi[0]);
    a.setYLim(m_dataLo[1], m_dataHi[1]);
    a.setXTicks({});
    a.setYTicks({});
    a.xaxis.ticks.gridOn = false;
    a.yaxis.ticks.gridOn = false;
    for (Spine& s : a.spines) s.visible = false;
    m_zmin = -radius;
    m_zmax = radius;
    setupDensity(fig, *tr);
    setRotation(fig, v.rotation);
    m_autoTitle = "Structure";

    try {
        m_shading = residueShading(ctx.project, ds, *tr);
    } catch (const std::out_of_range& e) {  // MDATK's KeyError("proximity")
        throw NoData(e.what());
    }
    m_resid = ints(required(*tr, "resid"));
    const double ps = v.point_size;
    const std::span<const Rgba> cycle = defaultColorCycle();
    const Rgba surface = toRgba(theme.surface);
    {
        Scatter s;
        s.zorder = 1;
        s.lineWidths = {0.0};
        s.faces = {cycle[0]};
        s.edges = s.faces;  // edgecolor "face"
        m_art.haloOuter = a.add(s);
        s.zorder = 2;
        s.faces = {cycle[1]};
        s.edges = s.faces;
        m_art.haloInner = a.add(s);
    }
    {
        Segments g;
        g.zorder = 3;
        g.capStyle = CapStyle::Round;
        g.colors = {cycle[0]};
        g.widths = {1.5};
        m_art.lines = a.add(g);
    }
    {
        Scatter s;
        s.zorder = 4;
        s.faces = {cycle[2]};
        s.edges = {surface};
        s.lineWidths = {0.6};
        m_art.points = a.add(s);
        Scatter ion;
        ion.zorder = 5;
        ion.sizes = {ps * 2.2};
        ion.faces = {toRgba(v.ion_color)};
        ion.edges = {surface};
        ion.lineWidths = {1.0};
        m_art.ions = a.add(ion);
    }
    if (v.show_time_label) {
        Text t;
        t.xy = {0.01, 0.99};
        t.coords = Coords::Axes;
        t.va = VAlign::Top;
        t.ha = HAlign::Left;
        t.color = toRgba(theme.ink2);
        t.font.size = ctx.fontSize * 0.85;
        m_art.timeText = a.add(t);
    }
    residueTags(fig);
    legend(fig);
    m_covered.assign(m_res, false);  // residues the 2D trace leaves to the 3D scene
    if (!v.representations.empty()) {
        setupStage(fig);
    }
}

void StructurePanel::viewCentre(const std::vector<Vec3d>& coords, std::size_t nFrames, std::size_t nRes) {
    const core::VisualParams& v = ctx.project.visual;
    m_shift.assign(nFrames, Vec3d{0.0, 0.0, 0.0});
    const auto modelCentre = [&] {
        // coords[0].mean(axis=0), of the finite points (RN-03)
        Vec3d sum{0.0, 0.0, 0.0};
        std::size_t n = 0;
        for (std::size_t r = 0; r < nRes; ++r) {
            if (!finite(coords[r])) continue;
            for (std::size_t c = 0; c < 3; ++c) sum[c] += coords[r][c];
            ++n;
        }
        if (n == 0) return Vec3d{0.0, 0.0, 0.0};
        for (double& c : sum) c /= static_cast<double>(n);
        return sum;
    };
    const std::string spec = trimmed(v.center_on);
    if (spec.empty()) {
        m_center = modelCentre();
        return;
    }
    if (lower(spec) == "origin") {
        m_center = {0.0, 0.0, 0.0};
        return;
    }
    std::size_t row = 0;
    try {
        const std::vector<std::int64_t> rows = residueRows(*ctx.dataset.result("trace"), spec);
        if (rows.empty()) {
            throw std::out_of_range("list index out of range");
        }
        row = static_cast<std::size_t>(rows[0]);
    } catch (const std::exception& exc) {  // SelectionError, IndexError
        m_viewNote = std::string("centre: ") + exc.what();
        warnOnce(std::string("structure view centre: ") + exc.what() + "; using the model centre");
        m_center = modelCentre();
        return;
    }
    if (v.center_follow) {
        for (std::size_t f = 0; f < nFrames; ++f) {
            for (std::size_t c = 0; c < 3; ++c) m_shift[f][c] = coords[f * nRes + row][c] - coords[row][c];
        }
    }
    m_center = coords[row];
}

std::vector<molecule::Vec3> StructurePanel::residueRgb(std::size_t pos) const {
    const core::VisualParams& v = ctx.project.visual;
    const std::vector<float> intensity = m_shading.colorIntensity.column(pos);
    const std::vector<Rgba> rgba = frameColors(m_shading.staticColors, v.color_mode, v.base_color, v.proximity_color,
                                               intensity, m_shading.member);
    std::vector<molecule::Vec3> out(rgba.size());
    const bool tint = !(m_stage && m_stage->renderer->drawsGlow());
    for (std::size_t k = 0; k < rgba.size(); ++k) {
        out[k] = {rgba[k].r, rgba[k].g, rgba[k].b};
        if (tint) {
            // numpy: the float32 glow column times Python floats stays float32
            const float g = std::clamp(m_shading.glow.at(k, pos) * static_cast<float>(v.glow_strength), 0.0f, 1.0f);
            const auto keep = static_cast<double>(1.0f - 0.35f * g);
            for (std::size_t c = 0; c < 3; ++c) {
                out[k][c] = out[k][c] * keep + m_shading.glowRgb[c] * 0.35 * static_cast<double>(g);
            }
        }
    }
    return out;
}

molecule::Shading StructurePanel::stageShading(std::size_t pos) const {
    const core::VisualParams& v = ctx.project.visual;
    const Rgba surface = toRgba(theme.surface);
    molecule::Shading s;
    s.fog_rgb = {surface.r, surface.g, surface.b};
    s.fog_strength = v.depth_cue ? 0.45 : 0.0;
    std::vector<double> glow(m_res);
    for (std::size_t k = 0; k < m_res; ++k) {
        const float g = m_shading.glow.at(k, pos) * static_cast<float>(v.glow_strength);
        glow[k] = static_cast<double>(std::clamp(g, 0.0f, 1.0f));
    }
    s.glow = std::move(glow);
    s.glow_rgb = m_shading.glowRgb;
    return s;
}

void StructurePanel::setRotation(Figure& fig, std::span<const double> angles) {
    const Mat3 r = eulerToMatrix(angles);
    m_viewMatrix = r;
    m_viewXyz.resize(m_centred.size());
    for (std::size_t i = 0; i < m_centred.size(); ++i) m_viewXyz[i] = rotate(m_centred[i], r);
    m_ionViewXyz.resize(m_ionsCentred.size());
    for (std::size_t i = 0; i < m_ionsCentred.size(); ++i) m_ionViewXyz[i] = rotate(m_ionsCentred[i], r);
    updateDensity(fig);
}

void StructurePanel::setFast(bool on) { m_fast = on; }

// ---- density map shadow --------------------------------------------------------------------------

void StructurePanel::setupDensity(Figure& fig, const AnalysisResult& trace) {
    const core::Project& project = ctx.project;
    const core::DensityParams& d = project.density;
    if (!core::densityActive(project)) {
        return;
    }
    if (metaFlag(trace.meta, "schematic")) {
        m_densityNote = "density map: load the real reference frame to see it";
        return;
    }
    if (d.transform.size() != 12) {
        m_densityNote = "density map: not fitted yet";
        return;
    }
    try {
        m_map = io::loadMapCached(d.map_path);
        // the placement follows the molecule if this reference frame differs from the one the map
        // was placed on (other start frame, cleaned vs raw files)
        std::vector<core::Vec3> reference(m_res);
        for (std::size_t r = 0; r < m_res; ++r) {
            for (std::size_t c = 0; c < 3; ++c) reference[r][c] = m_centred[r][c] + m_center[c];
        }
        const molecule::AnchoredPlacement placed = molecule::anchored(d.transform, d.anchor, reference);
        m_toMap = placed.transform;
        if (!placed.note.empty()) {
            warnOnce(placed.note);
        }
    } catch (const std::exception& exc) {
        m_map.reset();
        m_toMap.reset();
        m_densityNote = std::string("density map: ") + exc.what();
        core::logWarning(std::string("density map not drawn: ") + exc.what());
        return;
    }
    m_level = finiteMapLevel(*m_map, d.level);  // RN-33
    const Rgba c = toRgba(d.color);
    m_densityRgb = {c.r, c.g, c.b};
    Image im;
    im.rows = 2;
    im.cols = 2;
    im.rgba.assign(16, 0);
    im.extent = {0.0, 1.0, 0.0, 1.0};
    im.originUpper = false;
    im.interpolation = "bilinear";
    im.zorder = d.front ? 3.5 : 0.5;
    m_art.density = ax(fig).add(im);
    ax(fig).setXLim(m_dataLo[0], m_dataHi[0]);
    ax(fig).setYLim(m_dataLo[1], m_dataHi[1]);
}

void StructurePanel::setDensityTransform(Figure& fig, std::span<const double> values) {
    if (!m_map) {
        return;
    }
    m_toMap = core::RigidTransform::fromList(values);
    updateDensity(fig);
}

void StructurePanel::updateDensity(Figure& fig) {
    if (!ok() || !m_art.density || !m_map || !m_toMap) {
        return;
    }
    Image& im = artist<Image>(fig, *m_art.density);
    // the visible window (afterLayout widens it to the panel's aspect) plus a margin
    const Range xl = ax(fig).xaxis.limits, yl = ax(fig).yaxis.limits;
    double drift = 0.0;
    for (const Vec3d& s : m_shift) drift = std::max(drift, norm3(s));  // a followed residue moves the map
    const double mx = 0.1 * (xl.hi - xl.lo) + drift, my = 0.1 * (yl.hi - yl.lo) + drift;
    const molecule::ShadowWindow window{xl.lo - mx, xl.hi + mx, yl.lo - my, yl.hi + my};
    const int n = m_fast ? 100 : 200;
    std::optional<molecule::Shadow> sh;
    try {
        sh = molecule::projectDensity(m_map, *m_toMap, m_viewMatrix, m_center, m_level, window, n, n);
    } catch (const std::exception& exc) {  // e.g. an absurd stored placement
        core::logWarning(std::string("density map not drawn: ") + exc.what());
        sh.reset();
    }
    if (!sh) {
        im.visible = false;
        return;
    }
    const core::DensityParams& d = ctx.project.density;
    im.rows = static_cast<int>(sh->rows);
    im.cols = static_cast<int>(sh->cols);
    im.rgba = rgbaBytes(molecule::shadowRgba(*sh, m_densityRgb, d.opacity));
    im.extent = sh->extent;
    m_densityExtent = sh->extent;
    im.visible = true;
}

bool StructurePanel::densityAnimated() const {
    if (!m_art.density) {
        return false;
    }
    const bool moving = std::any_of(m_shift.begin(), m_shift.end(),
                                    [](const Vec3d& s) { return s[0] != 0.0 || s[1] != 0.0 || s[2] != 0.0; });
    return ctx.project.density.front || moving;
}

void StructurePanel::followDensity(Figure& fig, std::size_t pos) {
    if (!m_densityExtent || !m_art.density) {
        return;
    }
    const bool moving = std::any_of(m_shift.begin(), m_shift.end(),
                                    [](const Vec3d& s) { return s[0] != 0.0 || s[1] != 0.0 || s[2] != 0.0; });
    if (!moving) {
        return;
    }
    const Vec3d d = rotate(m_shift[pos], m_viewMatrix);
    const std::array<double, 4>& b = *m_densityExtent;
    artist<Image>(fig, *m_art.density).extent = {b[0] - d[0], b[1] - d[0], b[2] - d[1], b[3] - d[1]};
}

std::optional<std::vector<double>> StructurePanel::nudgeDensity(Figure& fig, std::array<double, 2> shift,
                                                                const std::optional<Mat3>& viewRotation) {
    if (!hasDensity() || !m_toMap) {
        return std::nullopt;
    }
    const Mat3& v = m_viewMatrix;
    const Mat3 vt = transpose(v);
    const core::RigidTransform& t = *m_toMap;
    core::RigidTransform move;
    if (viewRotation) {
        if (!m_densityCentre) {
            m_densityCentre = densityCentroid(*m_map, m_level);
        }
        const Vec3d c = t.inverse().apply(*m_densityCentre);  // the map centre in the model frame
        const Mat3 q = multiply(multiply(vt, *viewRotation), v);
        const Vec3d qc = rotate(c, q);
        move.rotation = q;
        move.translation = {c[0] - qc[0], c[1] - qc[1], c[2] - qc[2]};
    } else {
        move.translation = rotate({shift[0], shift[1], 0.0}, vt);
    }
    m_toMap = t.compose(move.inverse());
    updateDensity(fig);
    return m_toMap->toList();
}

std::vector<ArtistRef> StructurePanel::rotationArtists() const {
    if (!m_art.density || densityAnimated()) {  // already animated
        return {};
    }
    return {ref(*m_art.density)};
}

// ---- 3D representations ---------------------------------------------------------------------------

void StructurePanel::setupStage(Figure& fig) {
    const core::Dataset& ds = ctx.dataset;
    const core::Project& project = ctx.project;
    const core::VisualParams& v = project.visual;
    const AnalysisResult* atoms = ds.result("atoms");
    if (atoms == nullptr) {
        throw NoData("atoms");
    }
    const AnalysisResult& tr = *ds.result("trace");
    const NumArray& xyzArray = required(*atoms, "coords");
    const std::size_t nAtoms = xyzArray.shape().size() > 1 ? xyzArray.shape()[1] : 0;
    const std::size_t frames = xyzArray.shape().empty() ? 0 : xyzArray.shape()[0];
    molecule::TopologyInput input;
    input.names = metaStrings(atoms->meta, "atom_names");
    input.elements = metaStrings(atoms->meta, "elements");
    input.resnames = metaStrings(atoms->meta, "resnames");
    for (const std::int64_t r : ints(required(*atoms, "row"))) input.row.push_back(static_cast<std::int32_t>(r));
    input.n_res = static_cast<std::int32_t>(m_res);
    {
        const std::vector<std::int64_t> b = ints(required(*atoms, "bonds"));
        for (std::size_t i = 0; i + 1 < b.size(); i += 2) {
            input.bonds.push_back({static_cast<std::int32_t>(b[i]), static_cast<std::int32_t>(b[i + 1])});
        }
    }
    std::vector<float> xyz32;
    if (xyzArray.dtype() == NumArray::DType::Float32) {
        const std::span<const float> s = xyzArray.values<float>();
        xyz32.assign(s.begin(), s.end());
    } else {
        for (const double d : xyzArray.toDouble()) xyz32.push_back(static_cast<float>(d));
    }
    input.ref_xyz.resize(nAtoms);
    for (std::size_t i = 0; i < nAtoms && frames > 0; ++i) {
        input.ref_xyz[i] = {xyz32[3 * i], xyz32[3 * i + 1], xyz32[3 * i + 2]};
    }
    input.chain = ints(required(tr, "chain"));
    if (const NumArray* ri = tr.array("resindex")) {
        input.resindex = ints(*ri);
    } else {
        input.resindex.resize(m_res);
        std::iota(input.resindex.begin(), input.resindex.end(), 0);
    }
    input.backbone_atom = v.backbone_atom;
    const std::vector<std::string> groupNames = core::groupNames(project);
    std::vector<std::vector<bool>> member;
    if (const NumArray* gm = tr.array("group_member")) {
        const std::vector<double> m = gm->toDouble();
        const std::size_t g = gm->shape().empty() ? 0 : gm->shape()[0];
        for (std::size_t k = 0; k < g; ++k) {
            std::vector<bool> row(m_res);
            for (std::size_t r = 0; r < m_res; ++r) row[r] = m[k * m_res + r] != 0.0;
            member.push_back(std::move(row));
        }
    } else {
        const std::vector<std::int64_t> gi = ints(required(tr, "group_index"));
        for (std::size_t k = 0; k < groupNames.size(); ++k) {
            std::vector<bool> row(m_res);
            for (std::size_t r = 0; r < m_res; ++r) row[r] = gi[r] == static_cast<std::int64_t>(k);
            member.push_back(std::move(row));
        }
    }
    std::unique_ptr<molecule::MoleculeScene> scene =
        molecule::buildScene(input, member, groupNames, v.representations, v.rep_rest);
    if (!scene) {
        return;
    }
    StageRendererChoice choice = m_renderers ? m_renderers->choose(v.renderer) : defaultStageRenderer(v.renderer);
    if (!choice.renderer) {
        choice = defaultStageRenderer("software");
    }
    if (!choice.note.empty()) {
        m_stageNote = choice.note;
        core::logWarning("3D representations: " + choice.note);
    }
    auto stage = std::make_unique<Stage>();
    stage->scene = std::move(scene);
    stage->renderer = choice.renderer;
    stage->nAtoms = nAtoms;
    // atoms centred on the panel centre in float32, as MDATK computes them
    const std::array<float, 3> c32{static_cast<float>(m_center[0]), static_cast<float>(m_center[1]),
                                   static_cast<float>(m_center[2])};
    stage->atoms.resize(frames * nAtoms);
    for (std::size_t f = 0; f < frames; ++f) {
        const std::size_t sf = std::min(f, m_shift.size() - 1);
        const std::array<float, 3> s32{static_cast<float>(m_shift[sf][0]), static_cast<float>(m_shift[sf][1]),
                                       static_cast<float>(m_shift[sf][2])};
        for (std::size_t i = 0; i < nAtoms; ++i) {
            const std::size_t k = (f * nAtoms + i) * 3;
            for (std::size_t c = 0; c < 3; ++c) {
                const float x = (xyz32[k + c] - c32[c]) - s32[c];
                stage->atoms[f * nAtoms + i][c] = static_cast<double>(x);
            }
        }
    }
    Image im;
    im.fillsAxes = true;
    im.zorder = 3;
    im.rows = 0;
    im.cols = 0;
    stage->image = ax(fig).add(im);
    m_art.stage = stage->image;
    m_stage = std::move(stage);
    m_covered.assign(m_res, true);  // the scene draws every residue (or hides it): no 2D trace
}

void StructurePanel::renderStage(Figure& fig) {
    Stage& st = *m_stage;
    const std::size_t pos = *m_pos;
    if (st.pos != pos || st.items.empty()) {  // rotations reuse the frame's geometry
        try {
            const std::vector<molecule::Vec3> rgb = residueRgb(pos);
            st.items = st.scene->frame(std::span<const molecule::Vec3>(st.atoms).subspan(pos * st.nAtoms, st.nAtoms), rgb);
            if (m_nIons > 0) {
                molecule::Spheres sph;
                const Rgba ionColor = toRgba(ctx.project.visual.ion_color);
                molecule::Colored ions;
                for (std::size_t i = 0; i < m_nIons; ++i) {
                    sph.xyz.push_back(m_ionsCentred[pos * m_nIons + i]);
                    sph.radius.push_back(1.2);
                    sph.owner.push_back(-1);
                    sph.atom.push_back(-1);
                    ions.rgb.push_back({ionColor.r, ionColor.g, ionColor.b});
                }
                ions.prim = std::move(sph);
                st.items.push_back(std::move(ions));
            }
        } catch (const std::exception& exc) {  // never break a video for the 3D layer
            core::logWarning(std::string("3D representation not drawn: ") + exc.what());
            st.items.clear();
        }
    }
    st.pos = pos;
    Image& im = artist<Image>(fig, st.image);
    im.rows = im.cols = 0;
    im.rgba.clear();
    if (!fig.frozen()) {
        return;  // rendered once the layout is final (afterLayout)
    }
    const Box bb = fig.axesBox(mainAxes());
    const auto w = static_cast<int>(std::lround(bb.width())), h = static_cast<int>(std::lround(bb.height()));
    if (w < 2 || h < 2) {
        return;
    }
    molecule::Viewport vp;
    for (std::size_t i = 0; i < 3; ++i) vp.rotation[i] = m_viewMatrix[i];
    vp.xlim = {ax(fig).xaxis.limits.lo, ax(fig).xaxis.limits.hi};
    vp.ylim = {ax(fig).yaxis.limits.lo, ax(fig).yaxis.limits.hi};
    vp.width = w;
    vp.height = h;
    vp.zrange = {m_zmin, m_zmax};
    try {
        molecule::RgbaImage img = st.renderer->render(st.items, vp, stageShading(pos),
                                                      m_fast ? molecule::RenderQuality::Fast : molecule::RenderQuality::High);
        st.failed = false;
        im.rows = static_cast<int>(img.height);
        im.cols = static_cast<int>(img.width);
        im.rgba = std::move(img.rgba);
    } catch (const std::exception& exc) {  // never break a video for the 3D layer
        if (!st.failed) {
            core::logWarning(std::string("3D representation not drawn (2D trace shown instead): ") + exc.what());
        }
        st.failed = true;
    }
}

// ---- tags, legend, layout -------------------------------------------------------------------------

void StructurePanel::residueTags(Figure& fig) {
    const core::VisualParams& v = ctx.project.visual;
    m_tagRows = tagRows(*ctx.dataset.result("trace"), v);
    for (const std::int64_t k : m_tagRows) {
        Text t;
        t.text = std::to_string(m_resid.at(static_cast<std::size_t>(k)));
        t.xy = {0.0, 0.0};
        t.ha = HAlign::Center;
        t.va = VAlign::Center;
        t.color = toRgba(theme.ink2);
        t.font.size = ctx.fontSize * 0.7;
        t.zorder = 6;
        t.clipOn = true;
        m_art.tags.push_back(ax(fig).add(t));
    }
}

void StructurePanel::legend(Figure& fig) {
    const core::Project& project = ctx.project;
    const core::VisualParams& v = project.visual;
    const AnalysisResult& tr = *ctx.dataset.result("trace");
    const double ms = 7;
    std::vector<LegendEntry> entries;
    const auto add = [&](std::string text, Rgba color, double size, std::optional<double> alpha = std::nullopt) {
        LegendEntry e;
        e.text = std::move(text);
        e.handle.style = LineStyle::None;
        e.handle.marker = Marker::Circle;
        e.handle.markerSize = size;
        e.handle.color = color;
        e.handle.markerFace = color;
        e.handle.markerEdge = color;
        e.handle.alpha = alpha;
        e.handle.label = e.text;
        entries.push_back(std::move(e));
    };
    if (v.color_mode == "chain") {
        const std::vector<std::int64_t> chain = ints(required(tr, "chain"));
        const std::vector<std::string> names = chainNames(tr);
        std::set<std::int64_t> chains(chain.begin(), chain.end());
        for (const std::int64_t c : chains) {
            const auto k = static_cast<std::size_t>(std::find(chain.begin(), chain.end(), c) - chain.begin());
            // RN-09: the chain's name as the trace (and the sequence panel) gives it
            const std::string name = k < names.size() ? names[k] : std::string();
            add("Chain " + (name.empty() ? std::to_string(c + 1) : name), m_shading.staticColors[k], ms);
        }
    } else {
        const std::vector<std::int64_t> gi = ints(required(tr, "group_index"));
        for (std::size_t k = 0; k < project.groups.size(); ++k) {
            if (std::find(gi.begin(), gi.end(), static_cast<std::int64_t>(k)) != gi.end()) {
                add(project.groups[k].name, toRgba(project.groups[k].color), ms);
            }
        }
        if (v.color_mode == "proximity") {
            add("\xe2\x89\xa4 " + formatG(v.proximity_cutoff) + " \xc3\x85 of " + v.proximity_set, toRgba(v.proximity_color), ms);
        }
    }
    if (v.glow_mode != "none") {
        add("Glow (" + v.glow_mode + ")", toRgba(v.glow_color), ms * 1.6, 0.5);
    }
    if (m_nIons > 0) {
        add(ionLegendLabel(metaStrings(tr.meta, "ion_labels")), toRgba(v.ion_color), ms);  // RN-08
    }
    if (!entries.empty()) {
        Legend L;
        L.loc = LegendLoc::LowerRight;
        L.frameOn = true;
        L.entries = std::move(entries);
        styleLegend(L, legendStyle(theme, ctx.fontSize));
        ax(fig).legend = std::move(L);
        m_hasLegend = true;
    }
}

void StructurePanel::afterLayout(Figure& fig, const FontSet& /*fonts*/) {
    // Equal x / y scale: grow the data window to the panel's pixel aspect ratio
    fitAspect(fig);
    updateDensity(fig);  // the shadow covers the final window
    if (ok() && m_stage && m_pos) {
        drawFrame(fig, *m_pos);  // the 3D layer at the final size (and the 2D trace if it fails)
    }
}

void StructurePanel::fitAspect(Figure& fig) {
    if (!ok()) {
        return;
    }
    const Box bb = fig.axesBox(mainAxes());
    if (bb.width() <= 0 || bb.height() <= 0) {
        return;
    }
    std::array<double, 2> span{m_dataHi[0] - m_dataLo[0], m_dataHi[1] - m_dataLo[1]};
    const std::array<double, 2> mid{0.5 * (m_dataLo[0] + m_dataHi[0]), 0.5 * (m_dataLo[1] + m_dataHi[1])};
    const double boxAspect = bb.height() / bb.width();
    if (span[1] / span[0] < boxAspect) {
        span[1] = span[0] * boxAspect;
    } else {
        span[0] = span[1] / boxAspect;
    }
    ax(fig).setXLim(mid[0] - span[0] / 2, mid[0] + span[0] / 2);
    ax(fig).setYLim(mid[1] - span[1] / 2, mid[1] + span[1] / 2);
}

std::vector<ArtistRef> StructurePanel::animatedArtists() const {
    if (!ok()) {
        return {};
    }
    std::vector<ArtistRef> arts{ref(m_art.haloOuter), ref(m_art.haloInner), ref(m_art.lines)};
    const bool front = ctx.project.density.front;
    if (densityAnimated() && !front) {
        arts.insert(arts.begin(), ref(*m_art.density));  // moving with a followed residue, behind everything
    }
    if (m_art.stage) {
        arts.push_back(ref(*m_art.stage));
    }
    if (m_art.density && front) {
        arts.push_back(ref(*m_art.density));  // over the moving layers in videos too (zorder 3.5)
    }
    arts.push_back(ref(m_art.points));
    arts.push_back(ref(m_art.ions));
    for (const std::size_t t : m_art.tags) arts.push_back(ref(t));
    if (m_art.timeText) {
        arts.push_back(ref(*m_art.timeText));
    }
    if (m_hasLegend) {  // keep the legend above moving residues
        arts.push_back({mainAxes(), 0, true});
    }
    return arts;
}

// ---- one frame ------------------------------------------------------------------------------------

void StructurePanel::drawFrame(Figure& fig, std::size_t pos) {
    const core::VisualParams& v = ctx.project.visual;
    m_pos = pos;
    const std::span<const Vec3d> xyz = view(pos);
    if (m_art.density) {
        followDensity(fig, pos);
    }
    if (m_stage) {
        renderStage(fig);
    }
    // the 3D scene draws the residues, ions and glow; if it fails, the 2D trace does
    const bool threeD = m_stage && !m_stage->failed;
    const double span = std::max(m_zmax - m_zmin, 1e-6);
    std::vector<double> depth(m_res, 1.0), z(m_res);
    std::vector<std::size_t> rows;  // finite residues (RN-03)
    for (std::size_t k = 0; k < m_res; ++k) {
        z[k] = xyz[k][2];
        if (v.depth_cue) depth[k] = (xyz[k][2] - m_zmin) / span;
        if (m_finite[pos * m_res + k]) rows.push_back(k);
    }
    std::vector<Rgba> colors = frameColors(m_shading.staticColors, v.color_mode, v.base_color, v.proximity_color,
                                           m_shading.colorIntensity.column(pos), m_shading.member);
    if (v.depth_cue) {  // fog toward the surface for distant residues
        const Rgba surf = toRgba(theme.surface);
        for (std::size_t k = 0; k < m_res; ++k) {
            const double fog = 0.45 * (1.0 - depth[k]);
            colors[k].r = colors[k].r * (1 - fog) + surf.r * fog;
            colors[k].g = colors[k].g * (1 - fog) + surf.g * fog;
            colors[k].b = colors[k].b * (1 - fog) + surf.b * fog;
        }
    }
    const double ps = v.point_size;
    std::vector<std::size_t> order = stableOrder(rows, z);
    if (threeD) {  // residues drawn in 3D leave the 2D trace
        std::erase_if(order, [&](std::size_t k) { return static_cast<bool>(m_covered[k]); });
    }
    {
        Scatter& pts = artist<Scatter>(fig, m_art.points);
        pts.offsets.clear();
        pts.sizes.clear();
        pts.faces.clear();
        for (const std::size_t k : order) {
            pts.offsets.push_back({xyz[k][0], xyz[k][1]});
            pts.sizes.push_back(ps * (0.55 + 0.9 * depth[k]));
            pts.faces.push_back(colors[k]);
        }
    }
    {
        std::vector<std::array<std::int64_t, 2>> links;
        for (const auto& l : m_links) {
            const auto i = static_cast<std::size_t>(l[0]), j = static_cast<std::size_t>(l[1]);
            if (!m_finite[pos * m_res + i] || !m_finite[pos * m_res + j]) continue;
            if (threeD && (m_covered[i] || m_covered[j])) continue;
            links.push_back(l);
        }
        Segments& seg = artist<Segments>(fig, m_art.lines);
        seg.segments.clear();
        if (!links.empty()) {
            std::vector<double> segDepth(links.size());
            std::vector<std::size_t> idx(links.size());
            for (std::size_t s = 0; s < links.size(); ++s) {
                segDepth[s] = 0.5 * (depth[static_cast<std::size_t>(links[s][0])] + depth[static_cast<std::size_t>(links[s][1])]);
                idx[s] = s;
            }
            const std::vector<std::size_t> so = stableOrder(idx, segDepth);
            seg.colors.clear();
            seg.widths.clear();
            for (const std::size_t s : so) {
                const auto i = static_cast<std::size_t>(links[s][0]), j = static_cast<std::size_t>(links[s][1]);
                seg.segments.push_back({{xyz[i][0], xyz[i][1]}, {xyz[j][0], xyz[j][1]}});
                seg.colors.push_back({0.5 * (colors[i].r + colors[j].r), 0.5 * (colors[i].g + colors[j].g),
                                      0.5 * (colors[i].b + colors[j].b), 0.5 * (colors[i].a + colors[j].a)});
                seg.widths.push_back(v.line_width * (0.5 + 0.9 * segDepth[s]));
            }
        }
    }
    {
        std::vector<std::size_t> sel;
        std::vector<float> intensity(m_res, 0.0f);
        if (!(threeD && m_stage->renderer->drawsGlow())) {  // else the renderer draws the glow around the 3D residues
            for (const std::size_t k : rows) {
                intensity[k] = m_shading.glow.at(k, pos) * static_cast<float>(v.glow_strength);
                if (intensity[k] > 0.01f) sel.push_back(k);
            }
        }
        Scatter& outer = artist<Scatter>(fig, m_art.haloOuter);
        Scatter& inner = artist<Scatter>(fig, m_art.haloInner);
        outer.offsets.clear();
        inner.offsets.clear();
        if (!sel.empty()) {
            outer.sizes.assign(sel.size(), ps * 9.0);
            inner.sizes.assign(sel.size(), ps * 3.8);
            outer.faces.clear();
            inner.faces.clear();
            for (const std::size_t k : sel) {
                outer.offsets.push_back({xyz[k][0], xyz[k][1]});
                inner.offsets.push_back({xyz[k][0], xyz[k][1]});
                const float ao = std::clamp(0.20f * intensity[k], 0.0f, 1.0f);
                const float ai = std::clamp(0.45f * intensity[k], 0.0f, 1.0f);
                outer.faces.push_back(withAlphaOf(m_shading.glowRgb, static_cast<double>(ao)));
                inner.faces.push_back(withAlphaOf(m_shading.glowRgb, static_cast<double>(ai)));
            }
            outer.edges = outer.faces;
            inner.edges = inner.faces;
        }
    }
    {
        Scatter& ionPts = artist<Scatter>(fig, m_art.ions);
        ionPts.offsets.clear();
        if (!threeD) {  // else the ions are spheres in the 3D scene
            for (std::size_t i = 0; i < m_nIons; ++i) {
                const Vec3d& p = m_ionViewXyz[pos * m_nIons + i];
                if (finite(p)) ionPts.offsets.push_back({p[0], p[1]});
            }
        }
    }
    if (!m_art.tags.empty()) {  // push each label outward, away from the projected centre
        double cx = 0.0, cy = 0.0;
        for (const std::size_t k : rows) {
            cx += xyz[k][0];
            cy += xyz[k][1];
        }
        if (!rows.empty()) {
            cx /= static_cast<double>(rows.size());
            cy /= static_cast<double>(rows.size());
        }
        const double push = 0.07 * m_radius / std::max(v.zoom, 0.05) + 1.5;
        for (std::size_t t = 0; t < m_art.tags.size(); ++t) {
            const auto k = static_cast<std::size_t>(m_tagRows[t]);
            Text& text = artist<Text>(fig, m_art.tags[t]);
            if (!m_finite[pos * m_res + k]) {
                text.visible = false;
                continue;
            }
            text.visible = true;
            const double ox = xyz[k][0] - cx, oy = xyz[k][1] - cy;
            const double nrm = std::sqrt(ox * ox + oy * oy);
            double ux = 0.0, uy = 1.0;
            if (nrm > 1e-6) {
                ux = ox / std::max(nrm, 1e-6);
                uy = oy / std::max(nrm, 1e-6);
            }
            text.xy = {xyz[k][0] + ux * push, xyz[k][1] + uy * push};
        }
    }
    if (m_art.timeText) {
        artist<Text>(fig, *m_art.timeText).text = ctx.formatTime(pos);
    }
}

// ---- registration ---------------------------------------------------------------------------------

std::span<const PanelRegistration> structurePanels() {
    static const PanelRegistration table[] = {{"structure", &Panel::make<StructurePanel>}};
    return table;
}

}  // namespace mdbridge::render::panels
