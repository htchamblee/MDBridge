#pragma once
// The structure view of MDATK's figure (mdatk/render/panels.py StructurePanel,
// mdatk/render/molecule/stage.py, renderers.py): the projected backbone trace (points, links,
// two glow halos, ions, residue tags, time label, legend) in a bare square data window with equal
// aspect, the density map's shadow (behind or in front, following a followed centre), and the 3D
// representations rendered into an image drawn 1:1 into the axes (the scene of the molecule
// library, the built-in software renderer or another StageRenderer), which falls back to the 2D
// trace for a frame whose 3D render fails. Qt-free.
//
// Differences from MDATK v0.0.4.1 (docs/DIFFERENCES.md, "Rendering: structure view"):
// - RN-03: non-finite trace / ion coordinates (a blown-up frame) are left out of the view's extent
//   (MDATK: the radius became NaN and Matplotlib raised "Axis limits cannot be NaN or Inf"); one
//   warning names the frames; those residues are not drawn in those frames.
// - RN-08: the ion legend entry names the ions' element(s) (Mg²⁺, Na⁺, K⁺, Ca²⁺, ...; "Ions" for
//   a mixture or unknown ions) where MDATK always wrote "Mg²⁺".
// - RN-09: the chain legend names the chains as the trace does ("Chain A"; the chain's number only
//   when it has no name), where MDATK wrote "Chain 1", "Chain 2" (owner decision D11).
// - RN-16: visual.renderer "vtk" (and "vulkan") means the GPU renderer; without one the built-in
//   renderer is used with the note "the GPU renderer (Vulkan) is not available; using the built-in
//   renderer" (MDATK: "VTK is not installed (pip install vtk); ...").
// - RN-33: the density map's contour level (mean + level x rms) is computed from the finite voxels,
//   so a map with NaN / infinite voxels still shows its shadow (MDATK: the level became NaN, and
//   with RN-30 no shadow was drawn).

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/rigid_transform.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/molecule/scene.hpp"
#include "mdbridge/render/molecule/software.hpp"
#include "mdbridge/render/orientation.hpp"
#include "mdbridge/render/panels/panel.hpp"
#include "mdbridge/render/shading.hpp"

namespace mdbridge::render::panels {

// A renderer of the 3D layer (MDATK's renderer objects: render(items, viewport, shading, quality)).
// Implementations must be usable from several threads at once.
class StageRenderer {
public:
    virtual ~StageRenderer();
    // Throws on failure: the panel then shows the 2D trace for that frame.
    virtual molecule::RgbaImage render(const std::vector<molecule::Colored>& items, const molecule::Viewport& viewport,
                                       const molecule::Shading& shading, molecule::RenderQuality quality) const = 0;
    // The renderer draws the glow tint and halo itself (the 2D halos stay empty, colours untinted).
    virtual bool drawsGlow() const { return true; }
    // Class name for records and notes ("SoftwareRenderer").
    virtual std::string name() const = 0;
};

// The built-in software renderer (render/molecule/software.hpp).
class SoftwareStageRenderer final : public StageRenderer {
public:
    explicit SoftwareStageRenderer(molecule::SoftwareRenderOptions options = {}) : m_renderer(options) {}
    molecule::RgbaImage render(const std::vector<molecule::Colored>& items, const molecule::Viewport& viewport,
                               const molecule::Shading& shading, molecule::RenderQuality quality) const override {
        return m_renderer.render(items, viewport, shading, quality);
    }
    std::string name() const override { return "SoftwareRenderer"; }

private:
    molecule::SoftwareRenderer m_renderer;
};

struct StageRendererChoice {
    std::shared_ptr<const StageRenderer> renderer;
    std::string note;  // why another renderer than the requested one is used (logged as a warning)
};

// visual.renderer -> renderer (MDATK's get_renderer). The GUI / frame loop install one that knows
// the Vulkan renderer (setStageRendererFactory); defaultStageRenderer() knows only the software
// renderer.
class StageRendererFactory {
public:
    virtual ~StageRendererFactory();
    virtual StageRendererChoice choose(std::string_view name) const = 0;
};
// "auto" / "software" (and anything else): the software renderer; "vtk" / "vulkan": the software
// renderer with the note of RN-16.
StageRendererChoice defaultStageRenderer(std::string_view name);
// The factory structure panels built from now on use (process-wide; nullptr: defaultStageRenderer).
// Thread-safe.
void setStageRendererFactory(std::shared_ptr<const StageRendererFactory> factory);
std::shared_ptr<const StageRendererFactory> stageRendererFactory();
inline constexpr std::string_view kNoGpuRendererNote =
    "the GPU renderer (Vulkan) is not available; using the built-in renderer";

// map_level from the finite voxels (RN-33): mean + sigmas * rms of the voxels that are finite;
// io::mapLevel's value when all are; NaN when none is.
double finiteMapLevel(const io::DensityMap& map, double sigmas);

// A density-map placement made by hand on the frame a structure panel shows (MDATK's
// _on_map_moved, centre_density): `d` gets the transform (12 numbers, model -> map), the anchor
// (`model`, the panel's modelPoints of that frame, rounded as project files store them), fit_cc 0
// and `fitInfo`. The window's map drags and *Centre on model* store their placements with it.
void storeDensityPlacement(core::DensityParams& d, std::span<const double> transform, std::span<const Vec3d> model,
                           std::string fitInfo);

// The legend label of the ions (RN-08) from the trace's meta["ion_labels"] ("MG25", "NA101",
// "K+3", ...): "Mg²⁺", "Na⁺", "K⁺", "Ca²⁺", "Zn²⁺", "Mn²⁺", "Cl⁻" for one element, "Ions" for
// several or unknown ones (and without labels).
std::string ionLegendLabel(const std::vector<std::string>& ionLabels);

// The structure family's registration (registry.cpp).
std::span<const PanelRegistration> structurePanels();

class StructurePanel final : public Panel {
public:
    StructurePanel(const core::PanelSpec& spec, Context& ctx);
    ~StructurePanel() override;

    void afterLayout(Figure& fig, const FontSet& fonts) override;
    std::vector<ArtistRef> animatedArtists() const override;

    // set_rotation: re-orient the model (degrees about X, Y, Z; R = Rz Ry Rx); the shadow follows.
    void setRotation(Figure& fig, std::span<const double> angles);
    // fast quality while the GUI drags (coarser shadow, half-resolution 3D layer); takes effect at
    // the next update.
    void setFast(bool on);
    bool fast() const { return m_fast; }
    // The renderer of the 3D layer while fast() (the GUI's GPU renderer while dragging); nullptr:
    // the panel's own renderer.
    void setFastRenderer(std::shared_ptr<const StageRenderer> renderer);
    // A new visual.zoom without rebuilding the figure (the GUI's mouse wheel): the data window and
    // the shadow; the artists move at the next update.
    void setZoom(Figure& fig, double zoom);
    // set_density_transform: a new map placement (12 numbers, model -> map). Throws
    // core::ValueError for unusable numbers (RigidTransform::fromList).
    void setDensityTransform(Figure& fig, std::span<const double> values);
    // nudge_density: move the drawn map by `shift` (Å along the screen's x / y) or rotate it about
    // its density centre by `viewRotation` (screen axes); the new placement, or nullopt without a
    // map.
    std::optional<std::vector<double>> nudgeDensity(Figure& fig, std::array<double, 2> shift = {0.0, 0.0},
                                                    const std::optional<Mat3>& viewRotation = std::nullopt);
    bool hasDensity() const { return ok() && m_map != nullptr && m_art.density.has_value(); }
    // Static artists that still change when the view rotates (the GUI redraws them).
    std::vector<ArtistRef> rotationArtists() const;
    // The 3D layer of frame `pos` as the stage renderer is given it, for the GUI's 3D view: the
    // scene's primitives (ions included), the view rotation, the square data window (xlim / ylim),
    // the depth range and the shading; viewport width / height are left to the caller. nullopt
    // without 3D representations, when the panel shows "No data" or for a position past the
    // frames.
    struct StageFrame {
        std::vector<molecule::Colored> items;
        molecule::Viewport viewport;
        molecule::Shading shading;
    };
    std::optional<StageFrame> stageFrame(std::size_t pos);

    // ---- state (tests, records)
    double radius() const { return m_radius; }
    std::array<double, 2> dataLo() const { return m_dataLo; }
    std::array<double, 2> dataHi() const { return m_dataHi; }
    const std::string& viewNote() const { return m_viewNote; }
    const std::string& stageNote() const { return m_stageNote; }
    const std::string& densityNote() const { return m_densityNote; }
    const Vec3d& center() const { return m_center; }
    const std::vector<Vec3d>& frameShift() const { return m_shift; }
    const Mat3& viewMatrix() const { return m_viewMatrix; }
    // the view coordinates of the residues / ions in frame `pos`
    std::span<const Vec3d> view(std::size_t pos) const;
    std::span<const Vec3d> ionView(std::size_t pos) const;
    // The residues' points of frame `pos` as the panel places them, before the view rotation
    // (MDATK's st.centred[pos] + st.center): the anchor of a placement made on the shown frame.
    std::vector<Vec3d> modelPoints(std::size_t pos) const;
    std::optional<core::RigidTransform> densityPlacement() const { return m_toMap; }
    bool hasStage() const { return m_stage != nullptr; }
    std::string stageRenderer() const;
    std::string stageQuality() const { return m_fast ? "fast" : "high"; }
    bool stageFailed() const;
    // artist indices in the axes (nullopt when the panel has no such artist)
    struct Artists {
        std::size_t haloOuter = 0, haloInner = 0, lines = 0, points = 0, ions = 0;
        std::optional<std::size_t> timeText, density, stage;
        std::vector<std::size_t> tags;
    };
    const Artists& artists() const { return m_art; }
    // The residue colours of a frame for the 3D renderer (tinted toward the glow colour when the
    // renderer does not draw the glow).
    std::vector<molecule::Vec3> residueRgb(std::size_t pos) const;
    molecule::Shading stageShading(std::size_t pos) const;

protected:
    void setup(Figure& fig) override;
    void drawFrame(Figure& fig, std::size_t pos) override;

private:
    struct Stage;
    void viewCentre(const std::vector<Vec3d>& coords, std::size_t nFrames, std::size_t nRes);
    void setupDensity(Figure& fig, const core::AnalysisResult& trace);
    void setupStage(Figure& fig);
    void residueTags(Figure& fig);
    void legend(Figure& fig);
    void fitAspect(Figure& fig);
    void dataWindow();
    void updateDensity(Figure& fig);
    bool densityAnimated() const;
    void followDensity(Figure& fig, std::size_t pos);
    void renderStage(Figure& fig);
    void stageItems(std::size_t pos);

    std::size_t m_frames = 0, m_res = 0, m_nIons = 0;
    std::vector<Vec3d> m_centred, m_ionsCentred;  // frames x residues / ions
    std::vector<Vec3d> m_viewXyz, m_ionViewXyz;
    std::vector<std::array<std::int64_t, 2>> m_links;
    std::vector<Vec3d> m_shift;
    Vec3d m_center{};
    Mat3 m_viewMatrix{};
    double m_radius = 1.0, m_zmin = -1.0, m_zmax = 1.0;
    std::array<double, 2> m_dataLo{}, m_dataHi{};
    std::string m_viewNote, m_stageNote, m_densityNote;
    Shading m_shading;
    std::vector<std::int64_t> m_tagRows;
    std::vector<std::int64_t> m_resid;
    std::vector<bool> m_covered;
    std::vector<bool> m_finite;  // frames x residues: the trace point is finite
    Artists m_art;
    bool m_fast = false;
    double m_zoom = 1.0;
    std::shared_ptr<const StageRenderer> m_fastRenderer;
    std::optional<std::size_t> m_pos;
    std::shared_ptr<const StageRendererFactory> m_renderers;
    bool m_hasLegend = false;

    template <class T>
    T& artist(Figure& fig, std::size_t index) {
        return std::get<T>(fig.axes[mainAxes()].artists[index]);
    }
    ArtistRef ref(std::size_t index) const { return {mainAxes(), index, false}; }

    // density
    std::shared_ptr<const io::DensityMap> m_map;
    std::optional<core::RigidTransform> m_toMap;
    double m_level = 0.0;
    std::array<double, 3> m_densityRgb{};
    std::optional<std::array<double, 4>> m_densityExtent;
    std::optional<Vec3d> m_densityCentre;

    std::unique_ptr<Stage> m_stage;
};

}  // namespace mdbridge::render::panels
