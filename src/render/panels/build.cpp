// build_figure (build.hpp): mdatk/render/video.py build_figure and figure_size (MDATK's own code, ported).

#include "mdbridge/render/panels/build.hpp"

#include <algorithm>
#include <variant>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/log.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/python.hpp"

namespace mdbridge::render::panels {

std::pair<int, int> figureSize(int width, int height) {
    return {std::max(2, width / 2 * 2), std::max(2, height / 2 * 2)};
}

void BuiltFigure::update(std::size_t pos, const FontSet& fonts) {
    for (auto& p : panels) p->update(figure, pos);
    figure.draw(fonts);
}

std::vector<ArtistRef> BuiltFigure::animatedArtists() const {
    std::vector<ArtistRef> out;
    for (const auto& p : panels) {
        if (!p->ok()) continue;
        for (const ArtistRef& r : p->animatedArtists()) out.push_back(r);
    }
    return out;
}

BuiltFigure buildFigure(const core::Dataset& dataset, const core::Project& project, const FontSet& fonts, int width,
                        int height, double dpi, const std::string& watermark) {
    const core::VideoParams& vp = project.output.video;
    const auto [w, h] = figureSize(width > 0 ? width : static_cast<int>(vp.width),
                                   height > 0 ? height : static_cast<int>(vp.height));
    BuiltFigure out;
    out.context = std::make_unique<Context>(dataset, project);
    Context& ctx = *out.context;
    out.figure = Figure(w, h, dpi > 0 ? dpi : static_cast<double>(vp.dpi));
    Figure& fig = out.figure;
    fig.facecolor = toRgba(ctx.theme.page);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const core::LayoutPreset* preset = core::findLayoutPreset(project.layout.preset);
    if (!preset) {
        std::string names;
        for (const core::LayoutPreset& p : core::layoutPresets()) names += (names.empty() ? "" : ", ") + std::string(p.name);
        throw core::ValueError("unknown layout preset " + core::pythonRepr(core::Value(project.layout.preset)) +
                               " (choose from " + names + ")");
    }
    fig.grid.rows = preset->rows;
    fig.grid.cols = preset->cols;
    fig.grid.heightRatios.assign(preset->heightRatios.begin(), preset->heightRatios.end());
    const std::vector<core::PanelSpec>& specs = project.layout.panels;
    if (specs.size() > preset->slots.size()) {
        core::logWarning("layout " + project.layout.preset + " has " + std::to_string(preset->slots.size()) +
                         " slots; ignoring " + std::to_string(specs.size() - preset->slots.size()) + " extra panel(s)");
    }
    for (std::size_t k = 0; k < std::min(specs.size(), preset->slots.size()); ++k) {
        const core::LayoutSlot& s = preset->slots[k];
        const std::size_t first = fig.axes.size();
        out.panels.push_back(createPanel(fig, {s.row0, s.row1, s.col0, s.col1}, specs[k], ctx));
        // the axes the panel made, named for the records (panel / role)
        const Panel& p = *out.panels.back();
        for (std::size_t a = first; a < fig.axes.size(); ++a) {
            fig.axes[a].panel = static_cast<int>(k);
            fig.axes[a].role = a == p.mainAxes()                           ? "main"
                               : p.twinAxes() && a == *p.twinAxes()         ? "twin"
                               : p.colorbarAxes() && a == *p.colorbarAxes() ? "colorbar"
                                                                            : "other";
        }
    }
    if (!project.visual.title.empty()) {
        fig.setSuptitle(project.visual.title, toRgba(ctx.theme.ink), project.visual.font_size * 1.3);
    }
    if (!watermark.empty()) {
        Text& t = fig.addText(0.997, 0.003, watermark);
        t.ha = HAlign::Right;
        t.va = VAlign::Bottom;
        t.color = toRgba(ctx.theme.muted);
        t.font.size = project.visual.font_size * 0.75;
    }
    if (dataset.frameCount()) {
        for (auto& p : out.panels) p->update(fig, 0);
    }
    fig.draw(fonts);
    fig.freeze();  // positions stay: no re-layout (or jitter) per frame
    for (auto& p : out.panels) p->afterLayout(fig, fonts);
    fig.draw(fonts);
    return out;
}

void markAnimated(Figure& fig, const std::vector<ArtistRef>& artists, bool animated) {
    for (const ArtistRef& r : artists) {
        Axes& ax = fig.axes.at(r.axes);
        if (r.legend) {
            if (ax.legend) ax.legend->animated = animated;
            continue;
        }
        std::visit([animated](auto& a) { a.animated = animated; }, ax.artists.at(r.index));
    }
}

}  // namespace mdbridge::render::panels
