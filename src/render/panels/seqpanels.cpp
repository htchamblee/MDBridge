// The sequence and secondary-structure panels (panels/seqpanels.hpp), ported from MDATK's
// mdatk/render/seqpanels.py (MDATK's own code). The arithmetic follows MDATK's order of operations
// (float64, the glow in float32 as numpy computes it), so cell sizes, font sizes, positions and
// colours are MDATK's to the bit.

#include "mdbridge/render/panels/seqpanels.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::render::panels {

namespace {

constexpr double kNoteRow = 0.8;  // height of the 'choose a smaller range' line, in cells
constexpr std::int64_t kMaxLetterRows = 2;
constexpr Rgba kInterColor{1.0, 0.82, 0.0, 0.55};  // transparent yellow: pairs between chains

// Python's truthiness of a plain value.
bool truthy(const core::Value& v) {
    switch (v.type()) {
        case core::Value::Type::Null: return false;
        case core::Value::Type::Bool: return *v.get<bool>();
        case core::Value::Type::Int: return *v.get<std::int64_t>() != 0;
        case core::Value::Type::Float: return *v.get<double>() != 0.0;
        case core::Value::Type::String: return !v.get<std::string>()->empty();
        case core::Value::Type::Array:
        case core::Value::Type::Table: return v.size() != 0;
        case core::Value::Type::DateTime: return true;
    }
    return true;
}

// str(options.get(key, fallback) or ""), or str(options.get(key, fallback)) without `orEmpty`.
std::string optionText(const core::Table& options, std::string_view key, std::string_view fallback, bool orEmpty) {
    for (const auto& [k, v] : options) {
        if (k == key) {
            if (orEmpty && !truthy(v)) return {};
            return core::pythonStr(v);
        }
    }
    return std::string(fallback);
}

const core::Value* option(const core::Table& options, std::string_view key) {
    for (const auto& [k, v] : options) {
        if (k == key) return &v;
    }
    return nullptr;
}

// int(options.get("per_row", 0) or 0); RN-40: a value int() cannot read or a negative one is a
// "No data" message (MDATK: ValueError / TypeError traceback, or a reversed axis).
std::int64_t perRowOption(const core::Table& options) {
    const core::Value* v = option(options, "per_row");
    if (v == nullptr || !truthy(*v)) return 0;
    std::optional<std::int64_t> n;
    if (const auto* b = v->get<bool>()) {
        n = *b ? 1 : 0;
    } else if (const auto* i = v->get<std::int64_t>()) {
        n = *i;
    } else if (const auto* d = v->get<double>()) {
        if (std::isfinite(*d) && std::fabs(*d) < 9.2e18) n = static_cast<std::int64_t>(std::trunc(*d));
    } else if (const auto* s = v->get<std::string>()) {
        n = core::tryPythonInt(*s);
    }
    if (!n || *n < 0) {
        throw NoData("option 'per_row' must be a whole number of residues per row (0: automatic), not " +
                     core::pythonRepr(*v));
    }
    return *n;
}

// Python's float floor division a // b (float_floor_div of CPython's floatobject.c).
double pyFloorDiv(double a, double b) {
    double mod = std::fmod(a, b);
    double div = (a - mod) / b;
    if (mod != 0.0) {
        if ((b < 0) != (mod < 0)) {
            div -= 1.0;
        }
    }
    double floordiv = 0.0;
    if (div != 0.0) {
        floordiv = std::floor(div);
        if (div - floordiv > 0.5) floordiv += 1.0;
    } else {
        floordiv = std::copysign(0.0, a / b);
    }
    return floordiv;
}

std::int64_t pyMod(std::int64_t a, std::int64_t b) {
    const std::int64_t r = a % b;
    return (r != 0 && ((r < 0) != (b < 0))) ? r + b : r;
}
std::int64_t pyDiv(std::int64_t a, std::int64_t b) {
    const std::int64_t q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

Text makeText(std::string s, Point xy, HAlign ha, VAlign va, double size, double zorder) {
    Text t;
    t.text = std::move(s);
    t.xy = xy;
    t.coords = Coords::Data;
    t.ha = ha;
    t.va = va;
    t.font.size = size;
    t.zorder = zorder;
    return t;
}

Scatter emptyScatter(double zorder, std::string_view cycleColor) {
    Scatter s;
    s.zorder = zorder;
    const Rgba c = toRgba(cycleColor);  // Axes.scatter without colours: the next colour of the cycle
    s.faces = {c};
    s.edges = {c};  // edgecolors "face"
    s.lineWidths = {0.0};
    return s;
}

template <class T>
T& artist(Axes& ax, std::size_t i) {
    return std::get<T>(ax.artists.at(i));
}

}  // namespace

Rgba sequenceTextColor(const Rgba& c) {
    return toRgba(0.299 * c.r + 0.587 * c.g + 0.114 * c.b > 0.55 ? "#111111" : "#ffffff");
}

// ---- ResiduePanel ---------------------------------------------------------------------------------

ResiduePanel::ResiduePanel(std::string type, const core::PanelSpec& spec, Context& context)
    : Panel(std::move(type), spec, context) {}

template <class F>
void ResiduePanel::guarded(F&& body) {
    try {
        body();
    } catch (const NoData& e) {
        m_message = e.what();
        throw;
    } catch (const core::SelectionError& e) {  // residue_rows: KeyError(str(exc))
        m_message = e.what();
        throw NoData(m_message);
    } catch (const std::out_of_range& e) {  // a missing result or array: KeyError(name)
        m_message = e.what();
        throw NoData(m_message);
    }
}

void ResiduePanel::setup(Figure& fig) {
    guarded([&] { setupResidues(fig); });
}

void ResiduePanel::loadResidues() {
    const core::Dataset& ds = ctx.dataset;
    const core::AnalysisResult* tr = ds.result("trace");
    if (tr == nullptr) throw NoData("trace");
    if (ds.frameCount() == 0) throw NoData("no frames");  // RN-02
    const std::string range = optionText(spec().options, "residues", "", true);
    m_rows = residueRows(*tr, range);  // SelectionError -> "No data"
    const std::vector<std::string> names = residueNames(*tr);
    const std::vector<std::string> chains = chainNames(*tr);
    const core::NumArray* resid = tr->array("resid");
    if (resid == nullptr) throw std::out_of_range("resid");
    const std::vector<double> rv = resid->toDouble();
    // no nucleotide without a backbone atom (the detergent LDA; MDATK 0.1.x)
    const std::optional<std::vector<bool>> flags = backboneFlags(*tr);
    m_letters.clear();
    m_chains.clear();
    m_resid.clear();
    for (std::int64_t r : m_rows) {
        const auto k = static_cast<std::size_t>(r);
        m_letters.push_back(oneLetter(names.at(k), flags ? std::optional<bool>(flags->at(k)) : std::nullopt));
        m_chains.push_back(chains.at(k));
        m_resid.push_back(static_cast<std::int64_t>(rv.at(k)));
    }
    m_shading = residueShading(ctx.project, ds, *tr);
    const std::vector<std::int64_t> tags = tagRows(*tr, ctx.project.visual);
    const std::set<std::int64_t> tagSet(tags.begin(), tags.end());
    m_tagged.clear();
    for (std::int64_t r : m_rows) m_tagged.push_back(tagSet.count(r) != 0);
    m_rangeText = core::pythonStrip(range);
}

void ResiduePanel::bareAxes(Figure& fig) const {
    Axes& a = ax(fig);
    a.setXTicks({});
    a.setYTicks({});
    a.xaxis.ticks.gridOn = false;
    a.yaxis.ticks.gridOn = false;
    for (Spine& s : a.spines) s.visible = false;
}

std::vector<Rgba> ResiduePanel::colors(std::size_t pos) const {
    const core::VisualParams& v = ctx.project.visual;
    const std::vector<Rgba> all = frameColors(m_shading.staticColors, v.color_mode, v.base_color, v.proximity_color,
                                              m_shading.colorIntensity.column(pos), m_shading.member);
    std::vector<Rgba> out;
    out.reserve(m_rows.size());
    for (std::int64_t r : m_rows) out.push_back(all.at(static_cast<std::size_t>(r)));
    return out;
}

void ResiduePanel::setHalos(Figure& fig, const std::vector<Point>& xy, std::size_t pos, double size) const {
    // _glow(pos): float32 glow x glow_strength (numpy: a float32 array times a Python float)
    const float strength = static_cast<float>(ctx.project.visual.glow_strength);
    std::vector<std::size_t> sel;
    std::vector<float> intensity;
    for (std::size_t k = 0; k < m_rows.size(); ++k) {
        const float g = m_shading.glow.at(static_cast<std::size_t>(m_rows[k]), pos) * strength;
        if (g > 0.01f) {
            sel.push_back(k);
            intensity.push_back(g);
        }
    }
    Axes& a = ax(fig);
    Scatter& outer = artist<Scatter>(a, m_haloOuter);
    Scatter& inner = artist<Scatter>(a, m_haloInner);
    if (sel.empty()) {  // offsets only: sizes and colours stay as they were
        outer.offsets.clear();
        inner.offsets.clear();
        return;
    }
    const auto clip = [](float x) { return std::clamp(x, 0.0f, 1.0f); };
    std::vector<Point> at;
    std::vector<Rgba> fo, fi;
    for (std::size_t n = 0; n < sel.size(); ++n) {
        at.push_back(xy[sel[n]]);
        const Rgba base{m_shading.glowRgb[0], m_shading.glowRgb[1], m_shading.glowRgb[2], 0.0};
        Rgba o = base, i = base;
        o.a = clip(0.22f * intensity[n]);
        i.a = clip(0.45f * intensity[n]);
        fo.push_back(o);
        fi.push_back(i);
    }
    outer.offsets = at;
    outer.sizes.assign(sel.size(), size * 4.0);
    outer.faces = fo;
    outer.edges = fo;
    inner.offsets = at;
    inner.sizes.assign(sel.size(), size * 2.0);
    inner.faces = fi;
    inner.edges = fi;
}

// ---- sequence -------------------------------------------------------------------------------------

SequencePanel::SequencePanel(const core::PanelSpec& spec, Context& context)
    : ResiduePanel("sequence", spec, context) {}

void SequencePanel::setupResidues(Figure& fig) {
    loadResidues();
    bareAxes(fig);
    m_autoTitle = "Sequence" + (m_rangeText.empty() ? std::string() : " (" + m_rangeText + ")");
    m_perRowOption = perRowOption(spec().options);  // RN-40
    const std::size_t n = m_rows.size();
    // cells: residues, with a gap (and the next chain's letter) where the chain changes
    m_cells.clear();
    const std::set<std::string> distinct(m_chains.begin(), m_chains.end());
    const bool multi = distinct.size() > 1;
    for (std::size_t k = 0; k < n; ++k) {
        if (multi && (k == 0 || m_chains[k] != m_chains[k - 1])) {
            if (k != 0) m_cells.push_back({kGap, ""});
            m_cells.push_back({kChainLabel, m_chains[k].empty() ? "?" : m_chains[k]});
        }
        m_cells.push_back({static_cast<std::int64_t>(k), ""});
    }
    Axes& a = ax(fig);
    const Theme& th = theme;
    const double fs = ctx.fontSize;
    m_haloOuter = a.add(emptyScatter(1.0, "C0"));
    m_haloInner = a.add(emptyScatter(2.0, "C1"));
    Scatter boxes = emptyScatter(3.0, "C2");
    boxes.marker = Marker::Square;
    boxes.offsets.assign(n, Point{0.0, 0.0});
    boxes.sizes = {10.0};
    m_boxes = a.add(boxes);
    m_texts.clear();
    for (const std::string& ch : m_letters) {
        Text t = makeText(ch, {0.0, 0.0}, HAlign::Center, VAlign::Center, fs, 4.0);
        t.font.weight = 700;
        t.font.family = {"monospace"};
        m_texts.push_back(a.add(t));
    }
    m_tagTexts.clear();
    for (std::size_t k = 0; k < n; ++k) {
        if (!m_tagged[k]) continue;
        Text t = makeText(std::to_string(m_resid[k]), {0.0, 0.0}, HAlign::Center, VAlign::Top, fs * 0.6, 4.0);
        t.color = toRgba(th.ink2);
        m_tagTexts.emplace_back(k, a.add(t));
    }
    m_chainTexts.clear();
    for (const Cell& c : m_cells) {
        if (c.index != kChainLabel) continue;
        Text t = makeText(c.label, {0.0, 0.0}, HAlign::Center, VAlign::Center, fs * 0.8, 4.0);
        t.color = toRgba(th.muted);
        t.font.style = "italic";
        m_chainTexts.push_back(a.add(t));
    }
    Text note = makeText("", {0.995, 0.02}, HAlign::Right, VAlign::Bottom, fs * 0.65, 5.0);
    note.coords = Coords::Axes;
    note.color = toRgba(th.muted);
    m_note = a.add(note);
    m_xy.assign(n, Point{0.0, 0.0});
    m_cellSize = 10.0;
    layout(fig, 800.0, 60.0, 100.0, false);
}

void SequencePanel::afterLayout(Figure& fig, const FontSet& /*fonts*/) {
    if (!ok()) return;
    const Box b = fig.axesBox(mainAxes());
    if (b.width() > 0 && b.height() > 0) {
        layout(fig, b.width(), b.height(), fig.dpi(), true);  // the real size: warn now
    }
}

std::pair<std::vector<std::int64_t>, std::int64_t> SequencePanel::slots(std::int64_t perRow) const {
    std::vector<std::int64_t> out;
    std::int64_t k = 0;
    for (const Cell& c : m_cells) {
        if (c.index == kGap && pyMod(k, perRow) == 0) {
            out.push_back(-perRow);  // dropped gap: placed off the axes
            continue;
        }
        if (c.index == kChainLabel && perRow > 1 && pyMod(k, perRow) == perRow - 1) {
            ++k;
        }
        out.push_back(k);
        ++k;
    }
    const auto rows = static_cast<std::int64_t>(std::ceil(static_cast<double>(k) / static_cast<double>(perRow)));
    return {out, std::max<std::int64_t>(1, rows)};
}

void SequencePanel::layout(Figure& fig, double widthPx, double heightPx, double dpi, bool warn) {
    Axes& a = ax(fig);
    const double fs = ctx.fontSize;
    const double tagSpace = m_tagTexts.empty() ? 0.0 : 0.75;
    const double fontPx = fs * dpi / 72.0;
    const double noteSize = artist<Text>(a, m_note).font.size;
    const double notePx = 1.3 * noteSize * dpi / 72.0;  // the warning line, in pixels
    double cell = fontPx * 1.7;
    std::int64_t perRow = 1, nRows = 1;
    std::vector<std::int64_t> slotOf;
    double noteH = 0.0;
    while (true) {
        perRow = m_perRowOption != 0 ? m_perRowOption
                                     : std::max<std::int64_t>(1, static_cast<std::int64_t>(pyFloorDiv(widthPx, cell)));
        auto [s, r] = slots(perRow);
        slotOf = std::move(s);
        nRows = r;
        noteH = nRows > kMaxLetterRows ? std::max(kNoteRow, notePx / cell) : 0.0;
        if ((static_cast<double>(nRows) * (1.0 + tagSpace) + noteH) * cell <= heightPx || cell <= 6.0) break;
        cell *= 0.92;
    }
    if (m_perRowOption != 0) cell = std::min(cell, widthPx / static_cast<double>(perRow));
    m_perRow = perRow;
    m_nRows = nRows;
    const double rowH = 1.0 + tagSpace;
    a.setXLim(-0.5, static_cast<double>(perRow) - 0.5);
    a.setYLim(static_cast<double>(nRows) * rowH - 0.5 + noteH, -0.5);  // the note gets its own line
    const double pts = cell * 72.0 / dpi;  // cell edge in points
    m_cellSize = (0.9 * pts) * (0.9 * pts);
    const double letter = std::min(fs, 0.62 * pts);
    std::size_t chainK = 0;
    for (std::size_t c = 0; c < m_cells.size(); ++c) {
        const std::int64_t slot = slotOf[c];
        const double x = static_cast<double>(pyMod(slot, perRow));
        const double y = static_cast<double>(pyDiv(slot, perRow)) * rowH;
        const std::int64_t idx = m_cells[c].index;
        if (idx >= 0) {
            const auto k = static_cast<std::size_t>(idx);
            m_xy[k] = {x, y};
            Text& t = artist<Text>(a, m_texts[k]);
            t.xy = {x, y};
            t.font.size = letter;
            for (const auto& [row, art] : m_tagTexts) {
                if (row != k) continue;
                Text& tag = artist<Text>(a, art);
                tag.xy = {x, y + 0.55};
                tag.font.size = std::max(letter * 0.6, 4.0);
            }
        } else if (idx == kChainLabel) {
            Text& t = artist<Text>(a, m_chainTexts[chainK++]);
            t.xy = {x, y};
            t.font.size = letter * 0.85;
        }
    }
    Scatter& boxes = artist<Scatter>(a, m_boxes);
    boxes.offsets = m_xy;
    boxes.sizes.assign(m_xy.size(), m_cellSize);
    Text& note = artist<Text>(a, m_note);
    if (nRows > kMaxLetterRows) {
        const std::string msg = std::to_string(m_rows.size()) + " residues need " + std::to_string(nRows) +
                                " rows: choose a smaller range (option 'residues', e.g. A:1-30)";
        note.text = msg;
        if (warn) warnOnce("sequence panel: " + msg);
    } else {
        note.text.clear();
    }
}

std::vector<ArtistRef> SequencePanel::animatedArtists() const {
    if (!ok()) return {};
    std::vector<ArtistRef> out{ref(m_haloOuter), ref(m_haloInner), ref(m_boxes)};
    for (std::size_t i : m_texts) out.push_back(ref(i));
    for (const auto& [row, art] : m_tagTexts) out.push_back(ref(art));
    for (std::size_t i : m_chainTexts) out.push_back(ref(i));
    out.push_back(ref(m_note));
    return out;
}

void SequencePanel::drawFrame(Figure& fig, std::size_t pos) {
    const std::vector<Rgba> cols = colors(pos);
    Axes& a = ax(fig);
    Scatter& boxes = artist<Scatter>(a, m_boxes);
    boxes.faces = cols;
    boxes.edges = cols;
    for (std::size_t k = 0; k < m_texts.size(); ++k) {
        artist<Text>(a, m_texts[k]).color = sequenceTextColor(cols[k]);
    }
    setHalos(fig, m_xy, pos, m_cellSize);
}

// ---- secondary structure --------------------------------------------------------------------------

SecondaryStructurePanel::SecondaryStructurePanel(const core::PanelSpec& spec, Context& context)
    : ResiduePanel("secondary", spec, context) {}

void SecondaryStructurePanel::setupResidues(Figure& fig) {
    loadResidues();
    bareAxes(fig);
    m_source = optionText(spec().options, "source", "dynamic", false);
    m_autoTitle = "Secondary structure" + (m_rangeText.empty() ? std::string() : " (" + m_rangeText + ")");
    const std::size_t n = m_rows.size();
    const double gap = 1.5;
    m_x.assign(n, 0.0);
    for (std::size_t k = 1; k < n; ++k) {
        m_x[k] = m_x[k - 1] + (m_chains[k] != m_chains[k - 1] ? gap : 1.0);
    }
    m_xy.clear();
    for (double x : m_x) m_xy.push_back({x, 0.0});
    m_span = std::max(n ? m_x.back() - m_x.front() : 1.0, 1.0);
    m_rowOf.clear();
    for (std::size_t k = 0; k < n; ++k) m_rowOf.emplace_back(m_rows[k], k);
    std::sort(m_rowOf.begin(), m_rowOf.end());
    m_fixed.reset();
    m_pairmap.reset();
    if (m_source == "dotbracket") {
        m_fixed = dotbracketPairs(optionText(spec().options, "dotbracket", "", false));
    } else if (m_source == "dynamic" || m_source == "first") {
        const core::Dataset& ds = ctx.dataset;
        const std::vector<std::string> keys = ds.pairmapKeys();
        if (std::find(keys.begin(), keys.end(), "basepairs/map") == keys.end()) throw NoData("basepairs/map");
        m_pairmap = ds.pairmap("basepairs/map");
        if (m_source == "first") m_fixed = framePairs(0);
    } else {  // RN-41 (MDATK: any other value silently meant "dynamic")
        throw NoData("unknown source " + core::pythonRepr(core::Value(m_source)) +
                     " (option 'source': dynamic, first or dotbracket)");
    }
    const Theme& th = theme;
    const double fs = ctx.fontSize;
    m_markerSize = ctx.project.visual.point_size * 1.3;
    Axes& a = ax(fig);
    Line base;
    base.x = n ? std::vector<double>{m_x.front() - 0.6, m_x.back() + 0.6} : std::vector<double>{0.0, 1.0};
    base.y = {0.0, 0.0};
    base.color = toRgba(th.grid);
    base.width = 1.0;
    base.zorder = 1.0;
    a.add(base);
    m_haloOuter = a.add(emptyScatter(2.0, "C0"));
    m_haloInner = a.add(emptyScatter(3.0, "C1"));
    Segments intra;
    intra.colors = {toRgba(th.ink2)};
    intra.widths = {1.4};
    intra.zorder = 4.0;
    intra.capStyle = CapStyle::Round;
    m_intra = a.add(intra);
    Segments inter = intra;
    inter.colors = {kInterColor};
    inter.widths = {2.0};
    m_inter = a.add(inter);
    Scatter points = emptyScatter(5.0, "C2");
    points.offsets = m_xy;
    points.sizes = {m_markerSize};
    points.edges = {toRgba(th.surface)};
    points.lineWidths = {0.6};
    m_points = a.add(points);
    const bool showLetters = n <= 90;
    m_labels.clear();  // drawn over the arcs and halos, so animated as well
    for (std::size_t k = 0; k < n; ++k) {
        if (showLetters) {
            Text t = makeText(m_letters[k], {m_x[k], -0.06}, HAlign::Center, VAlign::Top, fs * 0.7, 5.0);
            t.color = toRgba(th.ink2);
            t.font.family = {"monospace"};
            m_labels.push_back(a.add(t));
        }
        if (m_tagged[k]) {
            Text t = makeText(std::to_string(m_resid[k]), {m_x[k], -0.16}, HAlign::Center, VAlign::Top, fs * 0.55, 5.0);
            t.color = toRgba(th.muted);
            m_labels.push_back(a.add(t));
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        if (k == 0 || m_chains[k] != m_chains[k - 1]) {
            Text t = makeText(m_chains[k], {m_x[k] - 0.6, 0.04}, HAlign::Right, VAlign::Bottom, fs * 0.75, 5.0);
            t.font.style = "italic";
            t.color = toRgba(th.muted);
            m_labels.push_back(a.add(t));
        }
    }
    a.setXLim((n ? m_x.front() : 0.0) - 1.2, (n ? m_x.back() : 1.0) + 0.8);
    a.setYLim(-1.05, 1.05);
    m_lastKey.reset();
}

std::vector<SecondaryStructurePanel::Pair> SecondaryStructurePanel::framePairs(std::size_t pos) const {
    if (!m_pairmap) return {};
    const core::PairMap& pm = *m_pairmap;
    const std::size_t frames = pm.frameCount();
    std::vector<Pair> out;
    if (frames == 0) return out;
    pos = std::min(pos, frames - 1);
    const auto s = static_cast<std::size_t>(pm.offsets.at(pos)), e = static_cast<std::size_t>(pm.offsets.at(pos + 1));
    const auto find = [&](std::int64_t row) -> std::optional<std::int64_t> {
        const auto it = std::lower_bound(m_rowOf.begin(), m_rowOf.end(), std::pair<std::int64_t, std::size_t>{row, 0});
        if (it == m_rowOf.end() || it->first != row) return std::nullopt;
        // a repeated trace row maps to its last panel row (a dict comprehension); rows are unique
        return static_cast<std::int64_t>(it->second);
    };
    for (std::size_t q = s; q < e && q < pm.i.size(); ++q) {
        const auto a = find(pm.i[q]), b = find(pm.j[q]);
        if (a && b && *a != *b) out.push_back({std::min(*a, *b), std::max(*a, *b), pm.values.at(q)});
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<SecondaryStructurePanel::Pair> SecondaryStructurePanel::dotbracketPairs(const std::string& text) const {
    if (core::pythonStrip(text).empty()) throw NoData("no dot-bracket given (option 'dotbracket')");
    chem::DotBracket db;
    try {
        db = chem::dotbracketToPairs(text);
    } catch (const core::ValueError& e) {
        throw NoData(std::string("dot-bracket: ") + e.what());
    }
    if (db.n != static_cast<std::int64_t>(m_rows.size())) {
        throw NoData("the dot-bracket has " + std::to_string(db.n) + " residues, the panel shows " +
                     std::to_string(m_rows.size()));
    }
    std::vector<Pair> out;
    for (const auto& [i, j] : db.pairs) out.push_back({i, j, 1.0});
    return out;
}

void SecondaryStructurePanel::arcs(const std::vector<Pair>& pairs, std::vector<std::vector<Point>>& intra,
                                   std::vector<std::vector<Point>>& inter) const {
    // np.linspace(0, pi, 25): i * (pi / 24), the last one pi itself
    constexpr int kPoints = 25;
    std::array<double, kPoints> theta{};
    const double step = std::numbers::pi / (kPoints - 1);
    for (int i = 0; i < kPoints; ++i) theta[static_cast<std::size_t>(i)] = static_cast<double>(i) * step;
    theta[kPoints - 1] = std::numbers::pi;
    for (const Pair& p : pairs) {
        const double x0 = m_x[static_cast<std::size_t>(p.a)], x1 = m_x[static_cast<std::size_t>(p.b)];
        const double r = 0.5 * (x1 - x0);
        const double height = std::max(r / (0.5 * m_span), 0.06);
        const bool between = m_chains[static_cast<std::size_t>(p.a)] != m_chains[static_cast<std::size_t>(p.b)];
        std::vector<Point> arc;
        arc.reserve(kPoints);
        for (double t : theta) {
            const double xs = 0.5 * (x0 + x1) - r * std::cos(t);
            const double ys = height * std::sin(t) * 0.95;
            arc.push_back(between ? Point{xs, -ys - 0.02} : Point{xs, ys + 0.02});
        }
        (between ? inter : intra).push_back(std::move(arc));
    }
}

std::vector<ArtistRef> SecondaryStructurePanel::animatedArtists() const {
    if (!ok()) return {};
    std::vector<ArtistRef> out{ref(m_haloOuter), ref(m_haloInner), ref(m_intra), ref(m_inter), ref(m_points)};
    for (std::size_t i : m_labels) out.push_back(ref(i));
    return out;
}

void SecondaryStructurePanel::drawFrame(Figure& fig, std::size_t pos) {
    const std::vector<Pair> pairs = m_fixed ? *m_fixed : framePairs(pos);
    Axes& a = ax(fig);
    if (!m_lastKey || *m_lastKey != pairs) {  // pairs change only every n frames
        std::vector<std::vector<Point>> intra, inter;
        arcs(pairs, intra, inter);
        artist<Segments>(a, m_intra).segments = std::move(intra);
        artist<Segments>(a, m_inter).segments = std::move(inter);
        m_lastKey = pairs;
    }
    artist<Scatter>(a, m_points).faces = colors(pos);
    setHalos(fig, m_xy, pos, m_markerSize * 2.2);
}

// ---- registry -------------------------------------------------------------------------------------

namespace {

// Panel::make, then the frame of every "No data" panel (RN-42: also when setup had hidden the
// spines) and the view of Matplotlib's cla() (limits 0..1; MDBridge's model would autoscale the
// empty data limits of the message's axes to +-0.055) with the grid off (rcParams axes.grid)
template <class P>
std::unique_ptr<Panel> makeResiduePanel(Figure& fig, SubplotSpec slot, const core::PanelSpec& spec, Context& ctx) {
    std::unique_ptr<Panel> p = Panel::make<P>(fig, slot, spec, ctx);
    if (!p->ok()) {
        Axes& a = fig.axes[p->mainAxes()];
        const AxesStyle style = axesStyle(ctx.theme, true, ctx.fontSize);
        a.spines[kLeft].visible = style.left.visible;
        a.spines[kBottom].visible = style.bottom.visible;
        a.spines[kRight].visible = style.right.visible;
        a.spines[kTop].visible = style.top.visible;
        a.setXLim(0.0, 1.0);
        a.setYLim(0.0, 1.0);
        a.xaxis.ticks.gridOn = false;
        a.yaxis.ticks.gridOn = false;
    }
    return p;
}

constexpr PanelRegistration kResiduePanels[] = {
    {"sequence", &makeResiduePanel<SequencePanel>},
    {"secondary", &makeResiduePanel<SecondaryStructurePanel>},
};

}  // namespace

std::span<const PanelRegistration> residuePanels() { return kResiduePanels; }

}  // namespace mdbridge::render::panels
