#pragma once
// MDATK's build_figure (mdatk/render/video.py): the figure of a project's layout with its panels,
// laid out and frozen, and the frame loop's split into a static background and animated artists.
// Qt-free.

#include <memory>
#include <string>
#include <vector>

#include "mdbridge/render/panels/panel.hpp"

namespace mdbridge::render::panels {

// The even pixel size of a video frame (figure_size: rounded down to even, at least 2).
std::pair<int, int> figureSize(int width, int height);

struct BuiltFigure {
    std::unique_ptr<Context> context;  // the panels refer to it
    Figure figure{2, 2, 72.0};
    std::vector<std::unique_ptr<Panel>> panels;  // in slot order

    // update(pos) of every panel, then Figure::draw(fonts) (the layout stays frozen)
    void update(std::size_t pos, const FontSet& fonts);
    // the animated artists of the panels that are ok(), in panel order
    std::vector<ArtistRef> animatedArtists() const;
};

// build_figure(dataset, project, width, height, dpi, watermark): the layout preset's grid (MDATK's
// constrained layout pads), a panel per slot (extra panels: "layout X has N slots; ignoring M extra
// panel(s)" warning), the suptitle (visual.title, 1.3 x font size), the watermark (bottom right),
// update(0) when there are frames, draw, freeze, after_layout. width / height / dpi <= 0: the
// project's output.video values. Throws core::ValueError for an unknown preset or panel type.
BuiltFigure buildFigure(const core::Dataset& dataset, const core::Project& project, const FontSet& fonts,
                        int width = 0, int height = 0, double dpi = 0.0, const std::string& watermark = {});

// Sets the `animated` flag of the given artists (render_video's set_animated(True) before the
// background is drawn): DisplayLayer::Static then leaves them out and DisplayLayer::Animated draws
// only them.
void markAnimated(Figure& fig, const std::vector<ArtistRef>& artists, bool animated = true);

}  // namespace mdbridge::render::panels
