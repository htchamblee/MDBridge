#pragma once
// The panels of MDATK's multi-panel figure (mdatk/render/panels.py, seqpanels.py; MDATK's own code,
// ported) on MDBridge's figure model (figure.hpp). Qt-free: the video loop, the charts, the preview
// and the GUI share them.
//
// The shared panel interface (wave 2: every panel type implements it):
//
//   Context   mdatk.render.panels.RenderContext: dataset, project, theme, the time axis in
//             inputs.time_unit, its label, the font size, format_time(pos)
//   Panel     one panel in one grid slot. Its life follows MDATK's build_figure / render_video:
//               1. createPanel(type, fig, slot, spec, ctx): add_subplot + style_axes, then setup()
//                  (the panel's axes, artists, legend, colour bar, twin axes); a NoData thrown by
//                  setup (MDATK's KeyError) turns the panel into the centred "No data: ..." text;
//                  then the title (spec.title or the automatic one, left, ink)
//               2. update(fig, 0) for every panel, Figure::draw(fonts), Figure::freeze(),
//                  afterLayout(fig, fonts) (structure: equal aspect; sequence: re-flow)
//               3. per video frame: update(fig, pos) moves the animated artists (the panel changes
//                  its artists in place; artist indices never change after setup), then
//                  Figure::draw(fonts) and buildDisplayList(fig, fonts, DisplayLayer::Animated)
//                  over the static background (DisplayLayer::Static, drawn once after
//                  markAnimated(fig, panels), build.hpp)
//   ArtistRef  which artist of which axes (or its legend) an animated artist is: the frame loop
//             sets their `animated` flag, the tests name them as the figure records do
//             ("ax<k>/<n>", "ax<k>/legend")
//   registry  createPanel / registeredPanelTypes (registry.cpp); each panel family contributes a
//             table of PanelRegistration (plotPanels(): timeseries, dual_axis, heatmap,
//             contact_map; structure; sequence and secondary)
//
// Panels keep axes and artist indices, never pointers or references into the Figure (its vectors
// grow while panels are built, and figures are copied or moved).

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/text.hpp"
#include "mdbridge/render/theme.hpp"

namespace mdbridge::render::panels {

// MDATK's KeyError in a panel's setup: the panel shows "No data: <message>" and the run continues.
class NoData : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// mdatk.render.panels.RenderContext.
struct Context {
    Context(const core::Dataset& dataset, const core::Project& project);

    const core::Dataset& dataset;
    const core::Project& project;
    const Theme& theme;
    std::vector<double> time;  // dataset.time(inputs.time_unit)
    std::string timeUnit;      // inputs.time_unit
    std::string timeLabel;     // "Time (ns)", "Frame", ...
    double fontSize = 11.0;    // visual.font_size

    // "frame N" (unit "frame") or "t = 1.234 ns  (frame N)" ("us" written as µs).
    std::string formatTime(std::size_t pos) const;

    // Dataset lookups as MDATK's panels make them: a missing key throws NoData with MDATK's
    // KeyError message ("analysis 'x' has not been computed", ...).
    core::Series series(std::string_view key) const;
    core::Matrix matrix(std::string_view key) const;
    core::PairMap pairmap(std::string_view key) const;
};

// An artist of the figure: artists[index] of axes `axes`, or that axes' legend.
struct ArtistRef {
    std::size_t axes = 0;
    std::size_t index = 0;
    bool legend = false;
    bool operator==(const ArtistRef&) const = default;
    // "ax<k>/<n>" or "ax<k>/legend" (the figure records' ids)
    std::string id() const;
};

class Panel {
public:
    virtual ~Panel() = default;
    Panel(const Panel&) = delete;
    Panel& operator=(const Panel&) = delete;

    const std::string& type() const { return m_type; }
    const core::PanelSpec& spec() const { return m_spec; }
    // false once the panel shows "No data: ..."
    bool ok() const { return m_ok; }
    const std::string& autoTitle() const { return m_autoTitle; }
    // the panel's axes: main, twin (dual_axis), colour bar (heatmap, contact_map)
    std::size_t mainAxes() const { return m_axes; }
    std::optional<std::size_t> twinAxes() const { return m_twin; }
    std::optional<std::size_t> colorbarAxes() const { return m_colorbar; }

    // The artists that change between frames (MDATK's animated_artists(); empty when not ok()).
    virtual std::vector<ArtistRef> animatedArtists() const { return {}; }
    // update(pos): moves the artists to frame position pos (nothing when not ok()).
    void update(Figure& fig, std::size_t pos) {
        if (m_ok) drawFrame(fig, pos);
    }
    // after_layout(): once the layout is frozen (axes sizes are final).
    virtual void afterLayout(Figure& /*fig*/, const FontSet& /*fonts*/) {}

    // Panel construction: add_subplot(slot), style_axes, setup() (NoData: "No data" text), title.
    // Used by the registry's factories: `PanelType` derives from Panel and has a constructor
    // (const core::PanelSpec&, Context&).
    template <class PanelType>
    static std::unique_ptr<Panel> make(Figure& fig, SubplotSpec slot, const core::PanelSpec& spec, Context& ctx) {
        std::unique_ptr<Panel> p(new PanelType(spec, ctx));
        p->init(fig, slot);
        return p;
    }

protected:
    Panel(std::string type, const core::PanelSpec& spec, Context& ctx);

    // The panel's content; throw NoData for missing data (MDATK's KeyError).
    virtual void setup(Figure& fig) = 0;
    virtual void drawFrame(Figure& /*fig*/, std::size_t /*pos*/) {}

    // the panel options (spec.options) read as MDATK's panels read them: options.get(key, default)
    // with Python's int() / float() / bool() of what is there
    const core::Value* option(std::string_view key) const;
    double optionNumber(std::string_view key, double fallback) const;
    std::int64_t optionInt(std::string_view key, std::int64_t fallback) const;
    bool optionBool(std::string_view key, bool fallback) const;
    std::string optionString(std::string_view key, std::string fallback) const;

    // the main axes of fig
    Axes& ax(Figure& fig) const { return fig.axes[m_axes]; }
    // style_axes(ax, theme, grid, font_size)
    void style(Axes& ax, bool grid = true) const;

    Context& ctx;
    const Theme& theme;
    std::string m_autoTitle;
    std::optional<std::size_t> m_twin, m_colorbar;

private:
    void init(Figure& fig, SubplotSpec slot);
    void showMessage(Figure& fig, const std::string& text);

    std::string m_type;
    core::PanelSpec m_spec;
    std::size_t m_axes = 0;
    bool m_ok = true;
};

using PanelFactory = std::unique_ptr<Panel> (*)(Figure& fig, SubplotSpec slot, const core::PanelSpec& spec, Context& ctx);

struct PanelRegistration {
    std::string_view type;  // the panel type name ("timeseries", ...)
    PanelFactory make;
};

// Each panel family's registrations (registry.cpp concatenates them in MDATK's registry order:
// structure, timeseries, dual_axis, heatmap, contact_map, sequence, secondary).
std::span<const PanelRegistration> plotPanels();  // timeseries, dual_axis, heatmap, contact_map (plots.cpp)

// The panel types that have an implementation, in registry order.
std::vector<std::string_view> registeredPanelTypes();
// panel_class(type)(fig, slot, spec, ctx). Throws core::ValueError("unknown panel type 'x' (choose
// from structure, timeseries, ...)") for a type no family registers (MDATK's message, listing
// core::panelTypes()).
std::unique_ptr<Panel> createPanel(Figure& fig, SubplotSpec slot, const core::PanelSpec& spec, Context& ctx);

}  // namespace mdbridge::render::panels
