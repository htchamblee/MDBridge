#pragma once
// The sequence and secondary-structure panels of MDATK's figure (mdatk/render/seqpanels.py; MDATK's
// own code, ported) on the shared panel interface (panel.hpp; registered as residuePanels()). Both
// show the residues of the structure view with the structure panel's colours, proximity glow and
// residue-number tags (shading.hpp), so a residue lights up in all three at once:
//
//   sequence   one-letter cells (bold monospace letters on boxes coloured like the structure, text
//              black or white by luminance), chain labels (italic) in the gaps between chains, tags
//              under the letters, wrapped into rows; the largest cell that fits the panel is chosen
//              again once the layout is final (afterLayout), and more than two rows add the note
//              "N residues need R rows: choose a smaller range ..." and log it once as a warning.
//              Options: residues ("A:1-20, B:5-30"), per_row (residues per row, 0: automatic).
//   secondary  residues along a line (a gap of 1.5 between chains), base pairs as arcs: within a
//              chain above the line (ink2), between chains transparent yellow below it; letters for
//              up to 90 residues, tags, chain names. Options: source ("dynamic": the basepairs
//              analysis' pairs of each frame, "first": its first frame's, "dotbracket": the option
//              dotbracket), dotbracket, residues.
//
// Each panel has one axes (bare: no ticks, grid or spines); missing data ("trace", "proximity",
// "basepairs/map", a range matching nothing, an unusable dot-bracket) gives MDATK's
// "No data: <reason>" message instead (ok() false; message() says why; the run continues).
//
// Differences from MDATK (docs/DIFFERENCES.md): an unusable per_row (negative, not a number)
// gives a "No data" message (MDATK: an inverted, garbled panel, or a ValueError traceback that
// aborts the render; RN-40); an unknown source is a "No data" message (MDATK: silently "dynamic";
// RN-41); a secondary panel's "No data" message keeps the frame every other "No data" panel has
// (MDATK: no frame once setup had hidden the spines; RN-42); a dataset without frames gives
// "No data: no frames" (RN-02).
//
// Qt-free.

#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "mdbridge/render/panels/panel.hpp"
#include "mdbridge/render/shading.hpp"

namespace mdbridge::render::panels {

// The registrations of this family: sequence, secondary (registry.cpp).
std::span<const PanelRegistration> residuePanels();

// The shared part of both panels (MDATK's _ResiduePanel).
class ResiduePanel : public Panel {
public:
    // The "No data" reason when !ok() (without the "No data: " prefix).
    const std::string& message() const { return m_message; }

protected:
    ResiduePanel(std::string type, const core::PanelSpec& spec, Context& ctx);

    // setup(): setupResidues() of the panel, where MDATK's KeyErrors (NoData, a residue range
    // matching nothing, a missing result or array) become NoData, their reason kept in message()
    void setup(Figure& fig) final;
    virtual void setupResidues(Figure& fig) = 0;
    template <class F>
    void guarded(F&& body);
    // _residues(): the selected trace rows, their letters, chains, residue numbers, shading, tags
    void loadResidues();
    // _bare_axes(): no ticks, grid or spines
    void bareAxes(Figure& fig) const;
    // _colors(pos): the rows' colours of the frame (frame_colors of the shading)
    std::vector<Rgba> colors(std::size_t pos) const;
    // _set_halos(xy, pos, size): the two glow halos of the rows that glow in the frame
    void setHalos(Figure& fig, const std::vector<Point>& xy, std::size_t pos, double size) const;
    ArtistRef ref(std::size_t index) const { return {mainAxes(), index, false}; }

    std::string m_message;
    std::vector<std::int64_t> m_rows;
    std::vector<std::string> m_letters;
    std::vector<std::string> m_chains;
    std::vector<std::int64_t> m_resid;
    Shading m_shading;
    std::vector<bool> m_tagged;
    std::string m_rangeText;
    std::size_t m_haloOuter = 0, m_haloInner = 0;
};

class SequencePanel final : public ResiduePanel {
public:
    // made by createPanel (Panel::make)
    SequencePanel(const core::PanelSpec& spec, Context& ctx);

    void afterLayout(Figure& fig, const FontSet& fonts) override;
    std::vector<ArtistRef> animatedArtists() const override;

    // A cell of the sequence: a residue (index into the rows), a gap between chains (kGap) or the
    // next chain's label (kChainLabel, with its text).
    static constexpr std::int64_t kGap = -1;
    static constexpr std::int64_t kChainLabel = -2;
    struct Cell {
        std::int64_t index = 0;
        std::string label;
    };
    const std::vector<Cell>& cells() const { return m_cells; }
    // _slots(per_row): the grid slot of every cell (a chain label never ends a row, a gap never
    // starts one: a dropped gap gets slot -perRow) and the number of rows.
    std::pair<std::vector<std::int64_t>, std::int64_t> slots(std::int64_t perRow) const;

    std::int64_t perRow() const { return m_perRow; }
    std::int64_t rows() const { return m_nRows; }
    double cellSize() const { return m_cellSize; }  // marker area of a box, points^2
    // the artists (indices into the axes): boxes, letters (per row), tags (per tagged row, in row
    // order), chain labels, the note
    std::size_t boxes() const { return m_boxes; }
    const std::vector<std::size_t>& letters() const { return m_texts; }
    const std::vector<std::pair<std::size_t, std::size_t>>& tags() const { return m_tagTexts; }  // (row, artist)
    const std::vector<std::size_t>& chainLabels() const { return m_chainTexts; }
    std::size_t note() const { return m_note; }

private:
    void setupResidues(Figure& fig) override;
    void drawFrame(Figure& fig, std::size_t pos) override;
    // _layout(width_px, height_px, dpi, warn): the largest cell that fits all rows; every artist placed
    void layout(Figure& fig, double widthPx, double heightPx, double dpi, bool warn);

    std::vector<Cell> m_cells;
    std::int64_t m_perRowOption = 0;
    std::int64_t m_perRow = 1;
    std::int64_t m_nRows = 1;
    double m_cellSize = 10.0;
    std::vector<Point> m_xy;
    std::size_t m_boxes = 0;
    std::vector<std::size_t> m_texts;
    std::vector<std::pair<std::size_t, std::size_t>> m_tagTexts;
    std::vector<std::size_t> m_chainTexts;
    std::size_t m_note = 0;
};

class SecondaryStructurePanel final : public ResiduePanel {
public:
    // made by createPanel (Panel::make)
    SecondaryStructurePanel(const core::PanelSpec& spec, Context& ctx);

    std::vector<ArtistRef> animatedArtists() const override;

    // A base pair as the panel draws it: row indices a < b (into the panel's rows) and the map value.
    struct Pair {
        std::int64_t a = 0, b = 0;
        double value = 1.0;
        auto operator<=>(const Pair&) const = default;
    };
    const std::string& source() const { return m_source; }
    double span() const { return m_span; }
    const std::vector<double>& positions() const { return m_x; }  // along the line, per row
    // _frame_pairs(pos): the map's pairs of the frame between rows of the panel, sorted, unique
    std::vector<Pair> framePairs(std::size_t pos) const;
    // the artists (indices into the axes)
    std::size_t intra() const { return m_intra; }
    std::size_t inter() const { return m_inter; }
    std::size_t points() const { return m_points; }
    const std::vector<std::size_t>& labels() const { return m_labels; }

private:
    void setupResidues(Figure& fig) override;
    void drawFrame(Figure& fig, std::size_t pos) override;
    std::vector<Pair> dotbracketPairs(const std::string& text) const;
    // _arcs(pairs): 25-point half-ellipses, height proportional to the span; within a chain above
    // the line, between chains below it
    void arcs(const std::vector<Pair>& pairs, std::vector<std::vector<Point>>& intra,
              std::vector<std::vector<Point>>& inter) const;

    std::string m_source;
    std::vector<double> m_x;
    std::vector<Point> m_xy;
    double m_span = 1.0;
    double m_markerSize = 0.0;
    std::vector<std::pair<std::int64_t, std::size_t>> m_rowOf;  // trace row -> panel row, sorted
    std::optional<std::vector<Pair>> m_fixed;
    std::optional<core::PairMap> m_pairmap;
    std::optional<std::vector<Pair>> m_lastKey;
    std::size_t m_intra = 0, m_inter = 0, m_points = 0;
    std::vector<std::size_t> m_labels;
};

// The text colour of a sequence letter on a box of colour `c`: "#111111" on light boxes (luminance
// 0.299 r + 0.587 g + 0.114 b above 0.55), else white.
Rgba sequenceTextColor(const Rgba& c);

}  // namespace mdbridge::render::panels
