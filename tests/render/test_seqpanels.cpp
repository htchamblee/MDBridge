// The sequence and secondary-structure panels (src/render/core seqpanels.*): counterparts of MDATK's
// tests/test_seqpanels.py (the panels' own tests; the structure-view tests are the structure
// panel's) and the panels' contents against MDATK's figures (data/render/figures/<case>.json,
// tests/parity/make_figure_cases.py; owner decision D10: contents only):
//
//   every case with a sequence / secondary panel is built from its recorded input dataset and
//   project; with MDATK's frozen axes boxes (the layout is the framework's subject, compared
//   within 1.5 % below) the panels re-flow after the layout as MDATK's do, and their records
//   (dumpFigure) equal MDATK's at every recorded position: artists in order with their data,
//   sizes, colours, texts, fonts and limits (numbers within 1e-9 relative, text extents within
//   a pixel), the animated artists, the panel attributes (per_row, n_rows, cell_size, source,
//   span) and the warnings logged.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/example.hpp"
#include "mdbridge/render/figure.hpp"
#include "mdbridge/render/figure_dump.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/panel.hpp"
#include "mdbridge/render/panels/seqpanels.hpp"
#include "mdbridge/render/shading.hpp"
#include "mdbridge/render/theme.hpp"
#include "figure_records.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using core::Array;
using core::Table;
using core::Value;
using namespace mdbridge::render;
using namespace mdbridge::render::panels;
namespace fs = std::filesystem;

namespace {

const fs::path kRender = test::kRenderData;

const FontSet& fonts() { return test::figureFonts(); }

const Value& at(const Value& v, std::string_view key) {
    const Value* p = v.find(key);
    if (!p) throw std::runtime_error("missing key " + std::string(key));
    return *p;
}
const Array& list(const Value& v) { return *v.get<Array>(); }
std::string str(const Value& v) { return v.isString() ? *v.get<std::string>() : std::string(); }
double num(const Value& v) {
    if (const auto* d = v.get<double>()) return *d;
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    return std::nan("");
}
bool isNum(const Value& v) { return v.get<double>() || v.get<std::int64_t>(); }

// tests/data/render/inputs/<name> (figure_records.hpp), read once
const core::Dataset& dataset(const std::string& name) {
    static std::map<std::string, std::unique_ptr<core::Dataset>> cache;
    auto& slot = cache[name];
    if (!slot) slot = std::make_unique<core::Dataset>(test::loadRenderInput(name));
    return *slot;
}

bool isResiduePanel(const std::string& type) { return type == "sequence" || type == "secondary"; }

// The record's main axes of panel k.
const Value* recordAxes(const Value& figure, int panel) {
    for (const Value& a : list(at(figure, "axes"))) {
        if (static_cast<int>(num(at(a, "panel"))) == panel && str(at(a, "role")) == "main") return &a;
    }
    return nullptr;
}

// build_figure(dataset, project, w, h, dpi, watermark) (panels/build.cpp) with the panel types that
// are registered; the others (until their chunks land) as bare stand-ins: their slot, styled axes
// and recorded title. With `recordedBoxes`, the frozen axes boxes are MDATK's (the record's), else
// MDBridge's layout's.
struct Built {
    std::unique_ptr<Context> context;
    Figure fig{2, 2, 72.0};
    std::vector<std::unique_ptr<Panel>> panels;  // per panel (nullptr: a stand-in)
    std::vector<std::size_t> axesOf;
};

bool registered(const std::string& type) {
    const std::vector<std::string_view> types = registeredPanelTypes();
    return std::find(types.begin(), types.end(), type) != types.end();
}

Built buildFigure(const Value& rec, const core::Dataset& ds, const core::Project& p, bool recordedBoxes) {
    const auto [w, h] = figureSize(static_cast<int>(num(list(at(rec, "size"))[0])),
                                   static_cast<int>(num(list(at(rec, "size"))[1])));
    Built b;
    b.context = std::make_unique<Context>(ds, p);
    b.fig = Figure(w, h, num(at(rec, "dpi")));
    Figure& fig = b.fig;
    const Theme& th = b.context->theme;
    fig.facecolor = toRgba(th.page);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const core::LayoutPreset* preset = core::findLayoutPreset(p.layout.preset);
    if (!preset) throw std::runtime_error("unknown preset " + p.layout.preset);
    fig.grid.rows = preset->rows;
    fig.grid.cols = preset->cols;
    fig.grid.heightRatios.assign(preset->heightRatios.begin(), preset->heightRatios.end());
    const Value& first = at(list(at(rec, "frames"))[0], "figure");
    const std::size_t n = std::min(p.layout.panels.size(), preset->slots.size());
    for (std::size_t k = 0; k < n; ++k) {
        const core::LayoutSlot& s = preset->slots[k];
        const SubplotSpec slot{s.row0, s.row1, s.col0, s.col1};
        const core::PanelSpec& spec = p.layout.panels[k];
        std::unique_ptr<Panel> panel;
        std::size_t a = 0;
        if (registered(spec.type)) {
            panel = createPanel(fig, slot, spec, *b.context);
            a = panel->mainAxes();
        } else {
            a = fig.addSubplot(slot);
            styleAxes(fig.axes[a], axesStyle(th, true, p.visual.font_size));
            const Value* ra = recordAxes(first, static_cast<int>(k));
            if (ra) {
                if (const Value* left = at(*ra, "titles").find("left")) {
                    Text& t = fig.axes[a].setTitle(str(at(*left, "text")), HAlign::Left);
                    t.color = toRgba(th.ink);
                    t.font.size = p.visual.font_size;
                }
            }
            if (spec.type == "structure") {
                fig.axes[a].setXTicks({});
                fig.axes[a].setYTicks({});
                for (Spine& sp : fig.axes[a].spines) sp.visible = false;
            }
        }
        fig.axes[a].panel = static_cast<int>(k);
        fig.axes[a].role = "main";
        b.axesOf.push_back(a);
        b.panels.push_back(std::move(panel));
    }
    if (!p.visual.title.empty()) fig.setSuptitle(p.visual.title, toRgba(th.ink), p.visual.font_size * 1.3);
    if (const std::string wm = str(at(rec, "watermark")); !wm.empty()) {
        Text& t = fig.addText(0.997, 0.003, wm);
        t.ha = HAlign::Right;
        t.va = VAlign::Bottom;
        t.color = toRgba(th.muted);
        t.font.size = p.visual.font_size * 0.75;
    }
    if (ds.frameCount()) {
        for (auto& panel : b.panels) {
            if (panel) panel->update(fig, 0);
        }
    }
    fig.draw(fonts());
    fig.freeze();
    if (recordedBoxes) {
        for (std::size_t k = 0; k < b.axesOf.size(); ++k) {
            const Value* ra = recordAxes(first, static_cast<int>(k));
            if (!ra) continue;
            const Array& px = list(at(*ra, "position_px"));
            const Box f{num(px[0]) / w, num(px[1]) / h, num(px[2]) / w, num(px[3]) / h};
            fig.axes[b.axesOf[k]].position = f;
            fig.axes[b.axesOf[k]].active = f;
        }
    }
    for (auto& panel : b.panels) {
        if (panel) panel->afterLayout(fig, fonts());
    }
    return b;
}

std::vector<std::string> caseNames() {
    std::vector<std::string> out;
    const fs::path dir = kRender / "figures";
    if (!fs::is_directory(dir)) return out;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.path().extension() != ".json") continue;
        const Value rec = test::readJson(e.path());
        bool mine = false;
        if (const Value* panels = rec.find("panels")) {
            for (const Value& pn : list(*panels)) mine = mine || isResiduePanel(str(at(pn, "type")));
        }
        if (mine) out.push_back(e.path().stem().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::size_t artistIndex(const std::string& id) { return static_cast<std::size_t>(std::stoul(id.substr(id.find('/') + 1))); }

// Compares a record of MDBridge's with MDATK's: every field MDATK recorded (but the walker's own
// "linestyles"), numbers within 1e-9 relative, text extents within `px` pixels.
class Comparer {
public:
    explicit Comparer(double px) : m_px(px) {}
    void value(const Value& mine, const Value& theirs, const std::string& where, bool extent = false) {
        if (isNum(theirs) && isNum(mine)) {
            const double a = num(mine), b = num(theirs);
            const double tol = extent ? m_px : 1e-9 * std::max(std::fabs(b), 1e-3);
            if (extent) m_maxExtent = std::max(m_maxExtent, std::fabs(a - b));
            if (!(std::fabs(a - b) <= tol) && !(std::isnan(a) && std::isnan(b))) fail(where, a, b);
            return;
        }
        if (theirs.isArray() && mine.isArray()) {
            const Array& ta = list(theirs);
            const Array& ma = list(mine);
            if (ta.size() != ma.size()) {
                m_failures.push_back(where + ": " + std::to_string(ma.size()) + " items, MDATK " + std::to_string(ta.size()));
                return;
            }
            for (std::size_t i = 0; i < ta.size(); ++i) value(ma[i], ta[i], where + "[" + std::to_string(i) + "]", extent);
            return;
        }
        if (theirs.isTable() && mine.isTable()) {
            for (const auto& [k, v] : *theirs.get<Table>()) {
                if (k == "linestyles" || k == "id") continue;
                const Value* m = mine.find(k);
                if (!m) {
                    m_failures.push_back(where + "." + k + ": missing");
                    continue;
                }
                value(*m, v, where + "." + k, k == "extent_px");
            }
            return;
        }
        if (!(mine == theirs)) {
            m_failures.push_back(where + ": " + core::toJson(mine, -1) + " vs MDATK " + core::toJson(theirs, -1));
        }
    }
    const std::vector<std::string>& failures() const { return m_failures; }
    double maxExtent() const { return m_maxExtent; }

private:
    void fail(const std::string& where, double a, double b) {
        char buf[96];
        std::snprintf(buf, sizeof buf, ": %.17g vs MDATK %.17g", a, b);
        m_failures.push_back(where + buf);
    }
    double m_px;
    double m_maxExtent = 0.0;
    std::vector<std::string> m_failures;
};

// a "No data" message is in axes coordinates; MDATK's walker names them "data" (transData equals
// transAxes for limits 0..1)
void normaliseCoords(Value& theirs, const Value& mine) {
    if (str(at(theirs, "kind")) != "text") return;
    if (str(at(theirs, "coords")) == "data" && str(at(mine, "coords")) == "axes" &&
        str(at(theirs, "text")).rfind("No data: ", 0) == 0) {
        *theirs.find("coords") = Value("axes");
    }
}

class SeqPanelCase : public ::testing::TestWithParam<std::string> {};

TEST_P(SeqPanelCase, ContentsAsMdatk) {
    const Value rec = test::readJson(kRender / "figures" / (GetParam() + ".json"));
    ASSERT_EQ(rec.find("error"), nullptr) << "MDATK failed in this case";
    const core::Project p = test::caseProject(rec);
    const std::string input = str(at(rec, "input"));
    const core::Dataset& ds = dataset(input);
    resetWarnOnce();
    test::LogCapture log;
    Built b = buildFigure(rec, ds, p, true);
    const Array& frames = list(at(rec, "frames"));
    const Value& first = at(frames[0], "figure");
    const Array& panels = list(at(rec, "panels"));
    // text extents: the plan's layout tolerance, 1.5 % of the figure size (the same fonts, but
    // Matplotlib's Agg measures hinted glyphs, MDBridge unhinted ones)
    Comparer c(0.015 * std::max(num(list(at(rec, "size"))[0]), num(list(at(rec, "size"))[1])));
    // the panels' own attributes
    for (std::size_t k = 0; k < b.panels.size(); ++k) {
        const auto* rp = dynamic_cast<const ResiduePanel*>(b.panels[k].get());
        if (!rp) continue;
        const ResiduePanel& pn = *rp;
        const Value& pr = panels[k];
        const std::string where = "panel " + std::to_string(k);
        EXPECT_EQ(pn.ok(), *at(pr, "ok").get<bool>()) << where << " " << pn.message();
        EXPECT_EQ(pn.autoTitle(), str(at(pr, "auto_title"))) << where;
        // the records name MDATK's axes; MDBridge's stand-ins have no twin / colour-bar axes
        std::vector<std::size_t> anim, mine;
        for (const Value& id : list(at(pr, "animated"))) anim.push_back(artistIndex(str(id)));
        for (const ArtistRef& r : pn.animatedArtists()) {
            EXPECT_EQ(r.axes, pn.mainAxes());
            EXPECT_FALSE(r.legend);
            mine.push_back(r.index);
        }
        EXPECT_EQ(mine, anim) << where;
        if (const auto* seq = dynamic_cast<const SequencePanel*>(&pn); seq && pn.ok()) {
            EXPECT_EQ(seq->perRow(), static_cast<std::int64_t>(num(at(pr, "per_row")))) << where;
            EXPECT_EQ(seq->rows(), static_cast<std::int64_t>(num(at(pr, "n_rows")))) << where;
            c.value(Value(seq->cellSize()), at(pr, "cell_size"), where + " cell_size");
        }
        if (const auto* sec = dynamic_cast<const SecondaryStructurePanel*>(&pn); sec && at(pr, "source").isString()) {
            if (pn.ok()) {
                EXPECT_EQ(sec->source(), str(at(pr, "source"))) << where;
                c.value(Value(sec->span()), at(pr, "span"), where + " span");
            }
        }
    }
    // the records at every position: the first in full, the later ones the animated artists' changes
    std::map<std::string, Value> firstRecords;  // artist id -> record at the first position
    for (const Value& ax : list(at(first, "axes"))) {
        for (const Value& a : list(at(ax, "artists"))) firstRecords[str(at(a, "id"))] = a;
    }
    for (std::size_t f = 0; f < frames.size(); ++f) {
        const std::size_t pos = static_cast<std::size_t>(num(at(frames[f], "pos")));
        // the fields that differ from the first position's
        std::map<std::string, Value> expected = firstRecords;
        if (f > 0) {
            for (const Value& change : list(at(frames[f], "artists"))) {
                Value& e = expected[str(at(change, "id"))];
                for (const auto& [k, v] : *change.get<Table>()) {
                    if (Value* field = e.find(k)) *field = v;
                }
            }
        }
        for (auto& panel : b.panels) {  // make_figure_cases: update(pos) at every position
            if (panel) panel->update(b.fig, pos);
        }
        b.fig.draw(fonts());
        const Value mine = dumpFigure(b.fig, fonts());
        for (std::size_t k = 0; k < b.panels.size(); ++k) {
            if (!dynamic_cast<const ResiduePanel*>(b.panels[k].get())) continue;
            const Value* ta = recordAxes(first, static_cast<int>(k));
            ASSERT_NE(ta, nullptr);
            const Value& ma = list(at(mine, "axes"))[b.axesOf[k]];
            const std::string where = "pos " + std::to_string(pos) + " panel " + std::to_string(k);
            if (f == 0) {
                for (const char* key : {"xlim", "ylim", "facecolor", "spines", "titles"}) {
                    c.value(at(ma, key), at(*ta, key), where + " " + key);
                }
                for (const char* axis : {"xaxis", "yaxis"}) {
                    c.value(at(at(ma, axis), "ticks"), at(at(*ta, axis), "ticks"), where + " " + axis + " ticks");
                    c.value(at(at(at(ma, axis), "tick_style"), "gridOn"), at(at(at(*ta, axis), "tick_style"), "gridOn"),
                            where + " " + axis + " grid");
                }
            }
            const Array& mArts = list(at(ma, "artists"));
            const Array& tArts = list(at(*ta, "artists"));
            ASSERT_EQ(mArts.size(), tArts.size()) << where << " artists";
            for (std::size_t n = 0; n < tArts.size(); ++n) {
                Value theirs = expected.at(str(at(tArts[n], "id")));
                normaliseCoords(theirs, mArts[n]);
                c.value(mArts[n], theirs, where + " artist " + std::to_string(n));
            }
        }
    }
    RecordProperty("max_text_extent_difference_px", std::to_string(c.maxExtent()));
    for (std::size_t i = 0; i < c.failures().size() && i < 30; ++i) ADD_FAILURE() << c.failures()[i];
    if (c.failures().size() > 30) ADD_FAILURE() << c.failures().size() << " differences in all";
    // the warnings of the panels (warn_once, once per case)
    std::vector<std::string> theirs, mine;
    for (const Value& l : list(at(rec, "log"))) {
        const std::string line = str(l);
        if (line.find("sequence panel:") != std::string::npos) theirs.push_back(line.substr(line.find(' ') + 1));
    }
    for (const std::string& l : log.lines()) {
        if (l.find("sequence panel:") != std::string::npos) mine.push_back(l.substr(l.find(' ') + 1));
    }
    EXPECT_EQ(mine, theirs);
}

INSTANTIATE_TEST_SUITE_P(Figures, SeqPanelCase, ::testing::ValuesIn(caseNames()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

// The panels' axes as MDBridge's constrained layout places them (with the setup-time layout of the
// sequence's letters, as MDATK lays the figure out before the re-flow), within 1.5 % of the figure
// size per edge, in the figures whose other panels are structure views (bare stand-ins here).
TEST(SeqPanels, LayoutAsMatplotlib) {
    int checked = 0;
    for (const std::string& name : caseNames()) {
        const Value rec = test::readJson(kRender / "figures" / (name + ".json"));
        const core::Project p = test::caseProject(rec);
        bool simple = true;
        for (const core::PanelSpec& s : p.layout.panels) simple = simple && (s.type == "structure" || isResiduePanel(s.type));
        if (!simple || !p.visual.title.empty()) continue;
        const std::string input = str(at(rec, "input"));
        Built b = buildFigure(rec, dataset(input), p, false);
        const Value& first = at(list(at(rec, "frames"))[0], "figure");
        const double tol = 0.015 * std::max(b.fig.widthPx(), b.fig.heightPx());
        for (std::size_t k = 0; k < b.panels.size(); ++k) {
            if (!dynamic_cast<const ResiduePanel*>(b.panels[k].get())) continue;
            const Value* ta = recordAxes(first, static_cast<int>(k));
            const Array& t = list(at(*ta, "position_px"));
            const Box m = b.fig.axesBox(b.axesOf[k]);
            EXPECT_NEAR(m.x0, num(t[0]), tol) << name << " panel " << k;
            EXPECT_NEAR(m.y0, num(t[1]), tol) << name << " panel " << k;
            EXPECT_NEAR(m.x1, num(t[2]), tol) << name << " panel " << k;
            EXPECT_NEAR(m.y1, num(t[3]), tol) << name << " panel " << k;
        }
        ++checked;
    }
    EXPECT_GE(checked, 10);
}

// ---- counterparts of MDATK's tests/test_seqpanels.py --------------------------------------------

core::Project projectWith(std::vector<core::PanelSpec> panels, const std::string& preset) {
    core::Project p;
    p.layout.preset = preset;
    p.layout.panels = std::move(panels);
    return p;
}

core::PanelSpec spec(const std::string& type, Table options = {}) { return {type, {}, "", std::move(options)}; }

// figure(project, size=(1600, 900)): build_figure of the example dataset at dpi 100 (the panel types
// not registered yet as bare stand-ins)
struct ExampleFigure {
    core::Dataset ds;
    core::Project p;
    std::unique_ptr<Context> context;
    std::unique_ptr<Figure> fig;
    std::vector<std::unique_ptr<Panel>> panels;
    std::vector<std::size_t> axesOf;
    template <class T>
    T& panel(std::size_t k) {
        return dynamic_cast<T&>(*panels.at(k));
    }
};

std::unique_ptr<ExampleFigure> figure(const core::Project& project, int w = 1600, int h = 900) {
    auto out = std::make_unique<ExampleFigure>(ExampleFigure{exampleDataset(project), project, nullptr, nullptr, {}, {}});
    out->context = std::make_unique<Context>(out->ds, out->p);
    out->fig = std::make_unique<Figure>(w, h, 100.0);
    Figure& fig = *out->fig;
    const Theme& th = out->context->theme;
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const core::LayoutPreset* preset = core::findLayoutPreset(project.layout.preset);
    fig.grid.rows = preset->rows;
    fig.grid.cols = preset->cols;
    fig.grid.heightRatios.assign(preset->heightRatios.begin(), preset->heightRatios.end());
    for (std::size_t k = 0; k < std::min(project.layout.panels.size(), preset->slots.size()); ++k) {
        const core::LayoutSlot& s = preset->slots[k];
        const SubplotSpec slot{s.row0, s.row1, s.col0, s.col1};
        const core::PanelSpec& spec = out->p.layout.panels[k];
        std::unique_ptr<Panel> panel = registered(spec.type) ? createPanel(fig, slot, spec, *out->context) : nullptr;
        const std::size_t a = panel ? panel->mainAxes() : fig.addSubplot(slot);
        if (!panel) {
            styleAxes(fig.axes[a], axesStyle(th, true, project.visual.font_size));
            fig.axes[a].setTitle("Stand-in", HAlign::Left);
        }
        out->axesOf.push_back(a);
        out->panels.push_back(std::move(panel));
    }
    for (auto& pn : out->panels) {
        if (pn) pn->update(fig, 0);
    }
    fig.draw(fonts());
    fig.freeze();
    for (auto& pn : out->panels) {
        if (pn) pn->afterLayout(fig, fonts());
    }
    fig.draw(fonts());
    return out;
}

template <class T>
const T& artistOf(const Figure& fig, const Panel& pn, std::size_t i) {
    return std::get<T>(fig.axes[pn.mainAxes()].artists.at(i));
}

TEST(SeqPanels, OneLetterCodes) {
    std::vector<std::string> got;
    for (const char* n : {"G5", "RA", "LG3", "0C", "U", "DT", "HIS", "XYZ"}) got.push_back(oneLetter(n));
    EXPECT_EQ(got, (std::vector<std::string>{"G", "A", "G", "C", "U", "T", "H", "X"}));
}

TEST(SeqPanels, ResidueRowsByChainAndNumber) {
    const core::Project p = projectWith({spec("sequence")}, "single");
    const core::Dataset ds = exampleDataset(p);
    const core::AnalysisResult& tr = *ds.result("trace");
    const std::vector<std::int64_t> rows = residueRows(tr, "B:13-15");
    const std::vector<double> resid = tr.array("resid")->toDouble();
    std::vector<double> got;
    for (std::int64_t r : rows) got.push_back(resid[static_cast<std::size_t>(r)]);
    EXPECT_EQ(got, (std::vector<double>{13, 14, 15}));
    EXPECT_EQ(residueRows(tr, "").size(), tr.array("chain")->size());
    EXPECT_THROW(residueRows(tr, "C:1"), core::SelectionError);
}

TEST(SeqPanels, StripLayoutsPutTheLastSlotAcrossTheBottom) {
    for (const core::LayoutPreset& preset : core::layoutPresets()) {
        const std::string name(preset.name);
        if (name.size() < 6 || name.substr(name.size() - 6) != "_strip") continue;
        ASSERT_EQ(preset.heightRatios.size(), static_cast<std::size_t>(preset.rows));
        const core::LayoutSlot& last = preset.slots.back();
        EXPECT_EQ(last.row0, preset.rows - 1);
        EXPECT_EQ(last.row1, preset.rows);
        EXPECT_EQ(last.col0, 0);
        EXPECT_EQ(last.col1, preset.cols);
        EXPECT_LT(preset.heightRatios.back(), 1.0);
        std::vector<core::PanelSpec> specs{spec("structure")};
        for (std::size_t k = 0; k + 2 < preset.slots.size(); ++k) specs.push_back({"timeseries", {"rmsd/Complex"}, "", {}});
        specs.push_back(spec("sequence"));
        auto f = figure(projectWith(specs, name));
        ASSERT_TRUE(f->panels.back()->ok()) << name;
        const Box strip = f->fig->axes[f->axesOf.back()].active, big = f->fig->axes[f->axesOf.front()].active;
        EXPECT_LT(strip.height(), big.height()) << name;
        EXPECT_GT(strip.width(), 0.8) << name;
    }
}

TEST(SeqPanels, SequenceWrapsAndWarns) {
    resetWarnOnce();
    test::LogCapture log;
    core::Project p = projectWith({spec("structure"), spec("sequence", Table{{"per_row", 5}})}, "single_strip");
    {
        auto f = figure(p);
        const auto& seq = f->panel<SequencePanel>(1);
        EXPECT_TRUE(seq.ok());
        EXPECT_GT(seq.rows(), 2);
        EXPECT_NE(artistOf<Text>(*f->fig, seq, seq.note()).text.find("choose a smaller range"), std::string::npos);
    }
    bool warned = false;
    for (const std::string& l : log.lines()) warned = warned || l.find("smaller range") != std::string::npos;
    EXPECT_TRUE(warned);
    p.layout.panels[1] = spec("sequence", Table{{"residues", "A:1-6"}});
    auto f = figure(p);
    auto& seq = f->panel<SequencePanel>(1);
    EXPECT_TRUE(seq.ok());
    EXPECT_EQ(seq.rows(), 1);
    EXPECT_EQ(artistOf<Text>(*f->fig, seq, seq.note()).text, "");
    ASSERT_EQ(seq.letters().size(), 6u);
    std::string letters;
    for (std::size_t i : seq.letters()) letters += artistOf<Text>(*f->fig, seq, i).text;
    EXPECT_EQ(letters, "GCAUGG");
    for (std::size_t pos : {0u, 5u}) seq.update(*f->fig, pos);  // colours and glow follow the frames
    EXPECT_EQ(artistOf<Scatter>(*f->fig, seq, seq.boxes()).faces.size(), 6u);
}

TEST(SeqPanels, SecondaryPanelPairsChangeWithTheFrames) {
    auto f = figure(projectWith({spec("secondary", Table{{"source", "dynamic"}})}, "single"));
    auto& arc = f->panel<SecondaryStructurePanel>(0);
    ASSERT_TRUE(arc.ok());
    arc.update(*f->fig, 0);
    const std::size_t n0 = artistOf<Segments>(*f->fig, arc, arc.inter()).segments.size();
    arc.update(*f->fig, 6);  // the example opens the duplex end every 6 frames
    const std::size_t n6 = artistOf<Segments>(*f->fig, arc, arc.inter()).segments.size();
    EXPECT_EQ(n0, 12u);
    EXPECT_EQ(n6, 10u);
    EXPECT_EQ(artistOf<Segments>(*f->fig, arc, arc.intra()).segments.size(), 0u);  // pairs between the two chains
    auto g = figure(projectWith({spec("secondary", Table{{"source", "first"}})}, "single"));
    auto& first = g->panel<SecondaryStructurePanel>(0);
    first.update(*g->fig, 6);
    EXPECT_EQ(artistOf<Segments>(*g->fig, first, first.inter()).segments.size(), 12u);
}

TEST(SeqPanels, SecondaryPanelFromDotbracket) {
    const std::string db = "((((....))))&" + std::string(12, '.');
    auto f = figure(projectWith({spec("secondary", Table{{"source", "dotbracket"}, {"dotbracket", db}})}, "single"));
    auto& arc = f->panel<SecondaryStructurePanel>(0);
    ASSERT_TRUE(arc.ok());
    arc.update(*f->fig, 0);
    EXPECT_EQ(artistOf<Segments>(*f->fig, arc, arc.intra()).segments.size(), 4u);
    EXPECT_EQ(artistOf<Segments>(*f->fig, arc, arc.inter()).segments.size(), 0u);
    auto bad = figure(projectWith({spec("secondary", Table{{"source", "dotbracket"}, {"dotbracket", "((..))"}})}, "single"));
    EXPECT_FALSE(bad->panels[0]->ok());  // 6 positions for 24 residues: a message instead of arcs
    EXPECT_EQ(bad->panel<ResiduePanel>(0).message(), "the dot-bracket has 6 residues, the panel shows 24");
}

TEST(SeqPanels, SecondaryPanelRequiresTheBasepairAnalysis) {
    // render.panels.required_analyses (analysis::panelRequiredAnalyses, milestone 3)
    core::Project p = projectWith({spec("secondary")}, "single");
    const auto has = [](const std::vector<std::string>& v, const char* name) {
        return std::find(v.begin(), v.end(), name) != v.end();
    };
    const std::vector<std::string> need = analysis::panelRequiredAnalyses(p);
    EXPECT_TRUE(has(need, "trace") && has(need, "basepairs"));
    p.layout.panels[0].options = Table{{"source", "dotbracket"}, {"dotbracket", "."}};
    EXPECT_FALSE(has(analysis::panelRequiredAnalyses(p), "basepairs"));
    p.layout.panels[0].options = Table{{"source", "dynamc"}};  // RN-41: still the analysis' (a message then)
    EXPECT_TRUE(has(analysis::panelRequiredAnalyses(p), "basepairs"));
    p.layout.panels = {spec("sequence")};
    p.visual.glow_mode = "none";
    p.visual.color_mode = "groups";
    EXPECT_EQ(analysis::panelRequiredAnalyses(p), (std::vector<std::string>{"trace"}));
    p.visual.color_mode = "proximity";
    EXPECT_EQ(analysis::panelRequiredAnalyses(p), (std::vector<std::string>{"proximity", "trace"}));
}

TEST(SeqPanels, LongChainNamesAreKept) {
    // the trace's chain names as they are (test_long_chain_names_are_kept: chain_names of the
    // structure; the panels show them)
    core::AnalysisResult tr;
    tr.name = "trace";
    tr.setArray("chain", core::NumArray(std::vector<std::int64_t>{0, 0, 1, 1}));
    tr.setArray("chain_name", core::NumArray(std::vector<std::string>{"RNA_chain_A", "RNA_chain_A", "RNA_chain_B", "RNA_chain_B"}));
    EXPECT_EQ(chainNames(tr), (std::vector<std::string>{"RNA_chain_A", "RNA_chain_A", "RNA_chain_B", "RNA_chain_B"}));
}

TEST(SeqPanels, LabelsDrawnOverAnimatedLayersAreAnimated) {
    core::Project p = projectWith({spec("structure"), spec("secondary", Table{{"source", "dynamic"}})}, "single_strip");
    auto f = figure(p);
    const auto& arc = f->panel<SecondaryStructurePanel>(1);
    const std::vector<ArtistRef> anim = arc.animatedArtists();
    ASSERT_FALSE(arc.labels().empty());
    for (std::size_t l : arc.labels()) {
        EXPECT_NE(std::find(anim.begin(), anim.end(), ArtistRef{arc.mainAxes(), l, false}), anim.end());
    }
    p.layout.panels[1] = spec("sequence");
    auto g = figure(p);
    const auto& seq = g->panel<SequencePanel>(1);
    const std::vector<ArtistRef> a2 = seq.animatedArtists();
    const auto in = [&](std::size_t i) {
        return std::find(a2.begin(), a2.end(), ArtistRef{seq.mainAxes(), i, false}) != a2.end();
    };
    for (std::size_t i : seq.chainLabels()) EXPECT_TRUE(in(i));
    for (const auto& [row, i] : seq.tags()) EXPECT_TRUE(in(i));
    EXPECT_TRUE(in(seq.note()));
    // and the frame loop draws them in the animated layer only
    markAnimated(*g->fig, a2);
    g->fig->draw(fonts());
    const DisplayList fg = buildDisplayList(*g->fig, fonts(), DisplayLayer::Animated);
    std::size_t texts = 0;
    for (const DisplayItem& it : fg.items) texts += std::holds_alternative<TextItem>(it.item) ? 1 : 0;
    EXPECT_GE(texts, seq.letters().size() + seq.chainLabels().size());
}

TEST(SeqPanels, SequenceChainLabelNeverEndsARow) {
    auto f = figure(projectWith({spec("structure"), spec("sequence")}, "single_strip"));
    const auto& seq = f->panel<SequencePanel>(1);
    for (std::int64_t perRow = 2; perRow < 20; ++perRow) {
        const auto [slots, rows] = seq.slots(perRow);
        ASSERT_EQ(slots.size(), seq.cells().size());
        for (std::size_t i = 0; i < slots.size(); ++i) {
            if (seq.cells()[i].index == SequencePanel::kChainLabel) {
                EXPECT_NE(slots[i] % perRow, perRow - 1);
            }
            if (seq.cells()[i].index == SequencePanel::kGap) {
                EXPECT_TRUE(slots[i] < 0 || slots[i] % perRow != 0);
            }
        }
        EXPECT_GE(rows * perRow, *std::max_element(slots.begin(), slots.end()) + 1);
    }
}

// ---- MDBridge's improvements (docs/DIFFERENCES.md) -------------------------------------------------

TEST(SeqPanels, UnusablePerRowIsAMessageRN40) {
    // MDATK: per_row = -3 inverts the x axis and gives boxes of 71 738 pt^2; "x" raises ValueError
    // out of build_figure (the whole render stops)
    for (const Value& v : {Value(std::int64_t{-3}), Value("x"), Value(-0.5 - 1.0)}) {
        auto f = figure(projectWith({spec("structure"), spec("sequence", Table{{"per_row", v}})}, "single_strip"));
        const ResiduePanel& seq = f->panel<ResiduePanel>(1);
        EXPECT_FALSE(seq.ok()) << core::toJson(v, -1);
        EXPECT_EQ(seq.message().rfind("option 'per_row' must be a whole number", 0), 0u) << seq.message();
        EXPECT_EQ(seq.autoTitle(), "Sequence");
    }
    // as MDATK: 0, empty, a fraction below 1 mean automatic; a numeric string or a float is truncated
    for (const Value& v : {Value(std::int64_t{0}), Value(""), Value(0.5), Value("7"), Value(7.9)}) {
        auto f = figure(projectWith({spec("structure"), spec("sequence", Table{{"per_row", v}})}, "single_strip"));
        const auto& seq = f->panel<SequencePanel>(1);
        EXPECT_TRUE(seq.ok()) << core::toJson(v, -1);
        if (num(v) == 7.9 || str(v) == "7") {
            EXPECT_EQ(seq.perRow(), 7);
        }
    }
}

TEST(SeqPanels, UnknownSourceIsAMessageRN41) {
    // MDATK: any source but "first" / "dotbracket" silently meant "dynamic"
    auto f = figure(projectWith({spec("secondary", Table{{"source", "dynamc"}})}, "single"));
    EXPECT_FALSE(f->panels[0]->ok());
    EXPECT_EQ(f->panel<ResiduePanel>(0).message(), "unknown source 'dynamc' (option 'source': dynamic, first or dotbracket)");
    EXPECT_EQ(f->panels[0]->autoTitle(), "Secondary structure");
}

TEST(SeqPanels, NoDataPanelsKeepTheirFrameRN42) {
    // MDATK: the secondary panel's message after its setup hid the spines had no frame
    auto f = figure(projectWith({spec("secondary", Table{{"source", "dotbracket"}})}, "single"));
    const Axes& ax = f->fig->axes[f->panels[0]->mainAxes()];
    EXPECT_FALSE(f->panels[0]->ok());
    EXPECT_TRUE(ax.spines[kLeft].visible);
    EXPECT_TRUE(ax.spines[kBottom].visible);
    EXPECT_FALSE(ax.spines[kTop].visible);
    EXPECT_FALSE(ax.spines[kRight].visible);
    EXPECT_EQ(ax.spines[kLeft].color, toRgba(lightTheme().grid));
}

TEST(SeqPanels, NoFramesIsAMessageRN02) {
    core::Project p = projectWith({spec("sequence"), spec("secondary", Table{{"source", "dotbracket"}, {"dotbracket", "."}})},
                                  "equal_2");
    const core::Dataset ex = exampleDataset(p);
    const core::Dataset empty({}, {}, ex.results(), ex.timeUnit());
    Figure fig(640, 360, 60.0);
    fig.grid.cols = 2;
    Context ctx(empty, p);
    for (int k = 0; k < 2; ++k) {
        auto pn = createPanel(fig, {0, 1, k, k + 1}, p.layout.panels[static_cast<std::size_t>(k)], ctx);
        EXPECT_FALSE(pn->ok());
        EXPECT_EQ(dynamic_cast<const ResiduePanel&>(*pn).message(), "no frames");
    }
}

TEST(SeqPanels, TextColourByLuminance) {
    EXPECT_EQ(sequenceTextColor(toRgba("#ffffff")), toRgba("#111111"));
    EXPECT_EQ(sequenceTextColor(toRgba("#2a78d6")), toRgba("#ffffff"));
    EXPECT_EQ(sequenceTextColor(toRgba("#ffb000")), toRgba("#111111"));
}

}  // namespace
