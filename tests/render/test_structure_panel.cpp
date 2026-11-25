// The structure view (src/render/panels structure.*, panel.*) against MDATK's figures
// (data/render/figures/<case>.json + .npz, tests/parity/make_figure_cases.py with MDBridge's
// improvements from tests/parity/improvements_render.py): for every case whose first panel is a
// structure view, the panel is built from the case's input dataset and project, the figure laid out,
// and the structure axes compared with the record at every recorded position - the 2D trace
// (points, links, halos, ions: offsets, sizes, colours, widths), tags, time label, legend, the
// density shadow (extent and RGBA) and the 3D layer (size and coverage) - plus the panel's own state
// (radius, window, notes, animated artists) and its warnings. Contents only (owner decision D10).
//
// The other panels of a case are wave 2's other chunks; here they are bare placeholder axes, and the
// structure axes is put at the box Matplotlib laid it out at (the framework's own layout is
// test_figure.cpp's subject), so the aspect fit, the density window and the 3D layer's size are
// compared exactly. Texts' and legends' extents are compared within the layout tolerance (1.5 % of
// the figure size).
//
// Also counterparts of MDATK's tests: test_seqpanels.py (zoom and centre, followed centre),
// test_density.py (shadow under the model, followed shadow, degraded placements),
// test_molecule_render.py (failed 3D render falls back to the 2D trace, front shadow animated above
// the 3D layer), and tests of the improvements RN-03, RN-08, RN-09, RN-16, RN-33.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/io/mrc.hpp"
#include "mdbridge/render/display_list.hpp"
#include "mdbridge/render/figure_dump.hpp"
#include "mdbridge/render/panels/build.hpp"
#include "mdbridge/render/panels/structure.hpp"
#include "mdbridge/render/shading.hpp"
#include "figure_records.hpp"

using namespace mdbridge;
using core::Array;
using core::Table;
using core::Value;
using namespace mdbridge::render;
using namespace mdbridge::render::panels;

namespace {

const FontSet& fonts() { return test::figureFonts(); }

// An array of a .npz archive written by render_common.write_npz (ZIP, members stored), widened to
// double: '<f4', '<f8', '<i4', '<i8', '|i1', '|u1', '|b1' in C order (analysis::readNpz has no
// uint8). Members by name without ".npy".
struct NpyArray {
    std::vector<std::size_t> shape;
    std::vector<double> values;
};

inline std::map<std::string, NpyArray> readStoredNpz(const std::filesystem::path& path) {
    const auto le32 = [](const unsigned char* p) {
        return std::uint32_t{p[0]} | (std::uint32_t{p[1]} << 8) | (std::uint32_t{p[2]} << 16) | (std::uint32_t{p[3]} << 24);
    };
    const auto le16 = [](const unsigned char* p) { return static_cast<std::uint16_t>(p[0] | (p[1] << 8)); };
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path.string());
    const std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::map<std::string, NpyArray> out;
    std::size_t pos = 0;
    while (pos + 30 <= bytes.size() && le32(bytes.data() + pos) == 0x04034b50) {
        const unsigned char* h = bytes.data() + pos;
        if (le16(h + 8) != 0) throw std::runtime_error(path.string() + ": compressed member");
        const std::size_t size = le32(h + 18), nl = le16(h + 26), el = le16(h + 28);
        std::string name(reinterpret_cast<const char*>(h) + 30, nl);
        const unsigned char* d = h + 30 + nl + el;
        if (name.size() > 4 && name.substr(name.size() - 4) == ".npy") name.resize(name.size() - 4);
        const std::size_t hl = d[6] == 1 ? le16(d + 8) : le32(d + 8);
        const std::size_t off = d[6] == 1 ? 10 : 12;
        const std::string header(reinterpret_cast<const char*>(d) + off, hl);
        const auto field = [&](const std::string& key) { return header.substr(header.find(':', header.find("'" + key + "'")) + 1); };
        std::string descr = field("descr");
        descr = descr.substr(descr.find('\'') + 1);
        descr = descr.substr(0, descr.find('\''));
        std::string sh = field("shape");
        sh = sh.substr(sh.find('(') + 1, sh.find(')') - sh.find('(') - 1);
        NpyArray a;
        std::size_t count = 1;
        std::stringstream ss(sh);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            if (tok.find_first_not_of(' ') != std::string::npos) {
                a.shape.push_back(static_cast<std::size_t>(std::stoull(tok)));
                count *= a.shape.back();
            }
        }
        const unsigned char* p = d + off + hl;
        a.values.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            if (descr == "<f4") {
                float f;
                std::memcpy(&f, p + 4 * i, 4);
                a.values[i] = f;
            } else if (descr == "<f8") {
                std::memcpy(&a.values[i], p + 8 * i, 8);
            } else if (descr == "<i4") {
                std::int32_t v;
                std::memcpy(&v, p + 4 * i, 4);
                a.values[i] = v;
            } else if (descr == "<i8") {
                std::int64_t v;
                std::memcpy(&v, p + 8 * i, 8);
                a.values[i] = static_cast<double>(v);
            } else if (descr == "|i1") {
                a.values[i] = static_cast<std::int8_t>(p[i]);
            } else if (descr == "|u1" || descr == "|b1") {
                a.values[i] = p[i];
            } else {
                throw std::runtime_error(path.string() + "/" + name + ": dtype " + descr);
            }
        }
        out.emplace(name, std::move(a));
        pos += 30 + nl + el + size;
    }
    return out;
}


const Value& at(const Value& v, std::string_view key) {
    const Value* f = v.find(key);
    if (f == nullptr) throw std::runtime_error("missing key " + std::string(key));
    return *f;
}
const Array& list(const Value& v) { return *v.get<Array>(); }
std::string str(const Value& v) { return v.stringOr(); }
double num(const Value& v) {
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    if (const auto* d = v.get<double>()) return *d;
    throw std::runtime_error("not a number");
}

// A panel of a type not registered in this build (wave 2's other chunks): a styled axes with its
// title.
class Placeholder final : public Panel {
public:
    Placeholder(const core::PanelSpec& spec, Context& ctx) : Panel(spec.type, spec, ctx) {}

protected:
    void setup(Figure& /*fig*/) override {}
};

// log capture: "LEVEL message" lines, as the generator's
struct LogCapture {
    std::vector<std::string> lines;
    std::mutex mutex;
    core::LogSink previous;
    LogCapture() {
        previous = core::setLogSink([this](core::LogLevel level, std::string_view m) {
            const std::lock_guard<std::mutex> lock(mutex);
            lines.push_back(std::string(core::logLevelName(level)) + " " + std::string(m));
        });
    }
    ~LogCapture() { core::setLogSink(previous); }
};

// build_figure(dataset, project, w, h, dpi, watermark) (build.hpp) with placeholders for the panel
// types not registered yet, the structure axes put at `structureBox` before after_layout
struct Built {
    std::unique_ptr<core::Dataset> dataset;
    std::unique_ptr<core::Project> project;
    std::unique_ptr<Context> ctx;
    std::unique_ptr<Figure> fig;
    std::vector<std::unique_ptr<Panel>> panels;
    StructurePanel* structure = nullptr;
    void update(std::size_t pos) {
        for (auto& p : panels) p->update(*fig, pos);
        fig->draw(fonts());
    }
};

Built buildFigure(core::Dataset dataset, core::Project project, int w, int h, double dpi,
                  const std::string& watermark = {}, const std::optional<Box>& structureBox = std::nullopt) {
    Built b;
    b.dataset = std::make_unique<core::Dataset>(std::move(dataset));
    b.project = std::make_unique<core::Project>(std::move(project));
    b.ctx = std::make_unique<Context>(*b.dataset, *b.project);
    std::tie(w, h) = figureSize(w, h);
    b.fig = std::make_unique<Figure>(w, h, dpi);
    Figure& fig = *b.fig;
    const Theme& th = b.ctx->theme;
    fig.facecolor = toRgba(th.page);
    fig.layout = {true, 0.08, 0.08, 0.06, 0.04};
    const core::LayoutPreset* preset = core::findLayoutPreset(b.project->layout.preset);
    if (preset == nullptr) throw std::runtime_error("unknown preset");
    fig.grid.rows = preset->rows;
    fig.grid.cols = preset->cols;
    fig.grid.heightRatios.assign(preset->heightRatios.begin(), preset->heightRatios.end());
    const std::vector<std::string_view> known = registeredPanelTypes();
    for (std::size_t k = 0; k < b.project->layout.panels.size() && k < preset->slots.size(); ++k) {
        const core::LayoutSlot& s = preset->slots[k];
        const core::PanelSpec& spec = b.project->layout.panels[k];
        const SubplotSpec slot{s.row0, s.row1, s.col0, s.col1};
        const std::size_t first = fig.axes.size();
        std::unique_ptr<Panel> p = std::find(known.begin(), known.end(), spec.type) != known.end()
                                       ? createPanel(fig, slot, spec, *b.ctx)
                                       : Panel::make<Placeholder>(fig, slot, spec, *b.ctx);
        for (std::size_t a = first; a < fig.axes.size(); ++a) {
            fig.axes[a].panel = static_cast<int>(k);
            fig.axes[a].role = a == p->mainAxes() ? "main" : "other";
        }
        if (auto* st = dynamic_cast<StructurePanel*>(p.get()); st && b.structure == nullptr) b.structure = st;
        b.panels.push_back(std::move(p));
    }
    if (!b.project->visual.title.empty()) fig.setSuptitle(b.project->visual.title, toRgba(th.ink), b.ctx->fontSize * 1.3);
    if (!watermark.empty()) {
        Text& t = fig.addText(0.997, 0.003, watermark);
        t.ha = HAlign::Right;
        t.va = VAlign::Bottom;
        t.color = toRgba(th.muted);
        t.font.size = b.project->visual.font_size * 0.75;
    }
    if (b.dataset->frameCount() > 0) {
        for (auto& p : b.panels) p->update(fig, 0);
    }
    fig.draw(fonts());
    fig.freeze();
    if (structureBox && b.structure != nullptr) {
        Axes& ax = fig.axes[b.structure->mainAxes()];
        const Box& q = *structureBox;
        ax.position = ax.active = {q.x0 / w, q.y0 / h, q.x1 / w, q.y1 / h};
        fig.draw(fonts());
    }
    for (auto& p : b.panels) p->afterLayout(fig, fonts());
    fig.draw(fonts());
    return b;
}

Box boxOf(const Value& v) {
    const std::vector<double> d = test::numbersOf(v);
    return {d.at(0), d.at(1), d.at(2), d.at(3)};
}

struct Compare {
    std::string where;
    double layoutTol = 0.0;  // pixels: extents
    int failures = 0;
    void fail(const std::string& what) {
        if (failures++ < 12) ADD_FAILURE() << where << ": " << what;
    }
    void value(const Value& mine, const Value& theirs, const std::string& path) {
        if (theirs.isNull() || mine.isNull()) {
            if (theirs.isNull() != mine.isNull()) fail(path + ": null vs not null");
            return;
        }
        if (const Table* t = theirs.get<Table>()) {
            if (!mine.isTable()) return fail(path + ": not a table");
            for (const auto& [k, v] : *t) {
                if (k == "linestyles" || k == "dtype" || k == "id") continue;
                const Value* m = mine.find(k);
                if (m == nullptr) {
                    fail(path + "." + k + ": missing");
                    continue;
                }
                if (k == "extent_px") {
                    const Box a = boxOf(*m), b = boxOf(v);
                    const double d = std::max({std::fabs(a.x0 - b.x0), std::fabs(a.x1 - b.x1), std::fabs(a.y0 - b.y0),
                                               std::fabs(a.y1 - b.y1)});
                    if (d > layoutTol) fail(path + ".extent_px off by " + std::to_string(d) + " px");
                    continue;
                }
                if (k == "coverage") {
                    const double a = num(*m), b = num(v);
                    if (std::fabs(a - b) > 0.005 * b + 20) fail(path + ".coverage " + std::to_string(a) + " vs " + std::to_string(b));
                    continue;
                }
                value(*m, v, path + "." + k);
            }
            return;
        }
        if (const Array* a = theirs.get<Array>()) {
            const Array* m = mine.get<Array>();
            if (m == nullptr) return fail(path + ": not an array");
            if (m->size() != a->size()) {
                return fail(path + ": " + std::to_string(m->size()) + " entries vs " + std::to_string(a->size()));
            }
            for (std::size_t i = 0; i < a->size(); ++i) value((*m)[i], (*a)[i], path + "[" + std::to_string(i) + "]");
            return;
        }
        if (theirs.isString() || mine.isString()) {
            if (str(mine) != str(theirs)) fail(path + ": '" + str(mine) + "' vs '" + str(theirs) + "'");
            return;
        }
        if (const bool* bt = theirs.get<bool>()) {
            const bool* bm = mine.get<bool>();
            if (bm == nullptr || *bm != *bt) fail(path + ": bool differs");
            return;
        }
        const double x = num(mine), y = num(theirs);
        if (std::isnan(x) && std::isnan(y)) return;
        if (!(std::fabs(x - y) <= 1e-9 * std::max(1.0, std::fabs(y)))) {
            char buf[160];
            std::snprintf(buf, sizeof buf, ": %.17g vs %.17g", x, y);
            fail(path + buf);
        }
    }
};

std::string refId(const ArtistRef& r) { return r.id(); }

std::vector<std::string> structureCases() {
    std::vector<std::string> out;
    const std::filesystem::path dir = test::kFigureCases;
    if (!std::filesystem::is_directory(dir)) return out;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".json") continue;
        const Value rec = test::untag(test::readJson(e.path()));
        const Value* panels = rec.find("project")->path({"layout", "panels"});
        if (panels && panels->isArray() && !list(*panels).empty() && str(at(list(*panels)[0], "type")) == "structure") {
            out.push_back(e.path().stem().string());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

class StructureCase : public ::testing::TestWithParam<std::string> {};

TEST_P(StructureCase, ContentsAsMdatk) {
    const Value rec = test::untag(test::readJson(test::kFigureCases / (GetParam() + ".json")));
    if (rec.find("error") != nullptr) {
        FAIL() << "MDATK raised for this case: " << str(at(at(rec, "error"), "message"))
               << " (MDBridge improves it: the expected data must not be an error)";
    }
    const std::vector<double> size = test::numbersOf(at(rec, "size"));
    const double dpi = num(at(rec, "dpi"));
    const Array& frames = list(at(rec, "frames"));
    ASSERT_FALSE(frames.empty());
    const Value& fig0 = at(at(frames[0], "figure"), "axes");
    const Value& panelRec = list(at(rec, "panels"))[0];
    const std::string axId = str(at(at(panelRec, "axes"), "main"));
    const Value* axRec = nullptr;
    for (const Value& a : list(fig0)) {
        if (str(at(a, "id")) == axId) axRec = &a;
    }
    ASSERT_NE(axRec, nullptr);

    resetWarnOnce();
    LogCapture log;
    Built b = buildFigure(test::loadRenderInput(str(at(rec, "input"))), test::caseProject(rec), static_cast<int>(size[0]),
                          static_cast<int>(size[1]), dpi, str(at(rec, "watermark")), boxOf(at(*axRec, "position_px")));
    ASSERT_NE(b.structure, nullptr);
    StructurePanel& st = *b.structure;
    const double tol = 0.015 * std::max(b.fig->widthPx(), b.fig->heightPx());

    // the panel record
    Compare cp{GetParam() + " panel", tol};
    EXPECT_EQ(st.ok(), *at(panelRec, "ok").get<bool>());
    EXPECT_EQ(st.autoTitle(), str(at(panelRec, "auto_title")));
    if (st.ok()) {
        Array animated;
        for (const ArtistRef& r : st.animatedArtists()) animated.push_back(refId(r));
        cp.value(Value(animated), at(panelRec, "animated"), "animated");
        cp.value(Value(st.radius()), at(panelRec, "radius"), "radius");
        cp.value(Value(Array{st.dataLo()[0], st.dataLo()[1]}), at(panelRec, "data_lo"), "data_lo");
        cp.value(Value(Array{st.dataHi()[0], st.dataHi()[1]}), at(panelRec, "data_hi"), "data_hi");
        EXPECT_EQ(st.viewNote(), str(at(panelRec, "view_note")));
        EXPECT_EQ(st.stageNote(), str(at(panelRec, "stage_note")));
        EXPECT_EQ(st.densityNote(), str(at(panelRec, "density_note")));
        if (const Value* stage = panelRec.find("stage")) {
            ASSERT_TRUE(st.hasStage());
            EXPECT_EQ(st.stageRenderer(), str(at(*stage, "renderer")));
            EXPECT_EQ(st.stageQuality(), str(at(*stage, "quality")));
            EXPECT_EQ(st.stageFailed(), *at(*stage, "failed").get<bool>());
        } else {
            EXPECT_FALSE(st.hasStage());
        }
    }

    // the images' pixels (density shadow) by npz member
    std::map<std::string, NpyArray> npz;
    const std::filesystem::path npzPath = test::kFigureCases / (GetParam() + ".npz");
    if (std::filesystem::exists(npzPath)) npz = readStoredNpz(npzPath);

    std::map<std::string, Value> first;  // pos 0's records by id, for later positions' changed fields
    for (const Value& fr : frames) {
        const auto pos = static_cast<std::size_t>(num(at(fr, "pos")));
        b.update(pos);
        const Value mine = dumpFigure(*b.fig, fonts());
        const Value* myAx = nullptr;
        for (const Value& a : list(at(mine, "axes"))) {
            if (str(at(a, "id")) == axId) myAx = &a;
        }
        ASSERT_NE(myAx, nullptr);
        std::map<std::string, const Value*> mineById;
        for (const Value& a : list(at(*myAx, "artists"))) mineById[str(at(a, "id"))] = &a;
        Compare c{GetParam() + " pos " + std::to_string(pos), tol};
        if (const Value* figure = fr.find("figure")) {
            (void)figure;
            c.value(at(*myAx, "xlim"), at(*axRec, "xlim"), "xlim");
            c.value(at(*myAx, "ylim"), at(*axRec, "ylim"), "ylim");
            c.value(at(*myAx, "position_px"), at(*axRec, "position_px"), "position_px");
            c.value(at(*myAx, "facecolor"), at(*axRec, "facecolor"), "facecolor");
            c.value(at(*myAx, "spines"), at(*axRec, "spines"), "spines");
            c.value(at(*myAx, "titles"), at(*axRec, "titles"), "titles");
            const Array& theirs = list(at(*axRec, "artists"));
            EXPECT_EQ(list(at(*myAx, "artists")).size(), theirs.size()) << GetParam() << ": artist count";
            for (const Value& a : theirs) {
                const std::string id = str(at(a, "id"));
                first[id] = a;
                if (!mineById.count(id)) {
                    c.fail(id + " missing");
                    continue;
                }
                c.value(*mineById[id], a, id);
            }
        } else {
            for (const Value& changed : list(at(fr, "artists"))) {
                const std::string id = str(at(changed, "id"));
                if (id.rfind(axId + "/", 0) != 0) continue;  // another panel's
                if (!mineById.count(id)) {
                    c.fail(id + " missing");
                    continue;
                }
                c.value(*mineById[id], changed, id);
            }
        }
        // images: the RGBA shown
        for (const auto& [id, rec2] : first) {
            const auto it = npz.find("p" + std::to_string(pos) + "/" + id + ".rgba");
            if (it == npz.end() || str(at(rec2, "kind")) != "image") continue;
            const std::size_t n = static_cast<std::size_t>(std::stoul(id.substr(id.find('/') + 1)));
            const Image& im = std::get<Image>(b.fig->axes[st.mainAxes()].artists.at(n));
            const std::vector<std::uint8_t> px = im.pixels();
            const std::vector<double>& want = it->second.values;
            ASSERT_EQ(px.size(), want.size()) << id;
            int bad = 0;
            for (std::size_t i = 0; i < px.size(); ++i) bad += std::fabs(px[i] - want[i]) > 1 ? 1 : 0;
            EXPECT_EQ(bad, 0) << c.where << " " << id << ": RGBA bytes off by more than 1";
        }
    }

    // warnings: the structure view's ones, as the record's log
    std::set<std::string> mineLog(log.lines.begin(), log.lines.end());
    for (const Value& l : list(at(rec, "log"))) {
        const std::string line = str(l);
        const bool ours = line.find("structure view") != std::string::npos || line.find("3D represent") != std::string::npos ||
                          line.find("density") != std::string::npos;
        if (ours) {
            EXPECT_TRUE(mineLog.count(line)) << GetParam() << ": expected log line " << line;
        }
    }
    std::set<std::string> theirLog;
    for (const Value& l : list(at(rec, "log"))) theirLog.insert(str(l));
    for (const std::string& l : mineLog) EXPECT_TRUE(theirLog.count(l)) << GetParam() << ": unexpected log line " << l;

    // the display list draws the 3D layer 1:1 into the axes box
    if (st.artists().stage) {
        const Image& im = std::get<Image>(b.fig->axes[st.mainAxes()].artists[*st.artists().stage]);
        const DisplayList dl = buildDisplayList(*b.fig, fonts());
        bool found = false;
        for (const DisplayItem& d : dl.items) {
            if (const auto* item = std::get_if<ImageItem>(&d.item); item && item->width == im.cols && item->height == im.rows) {
                const Box bb = b.fig->axesBox(st.mainAxes());
                EXPECT_DOUBLE_EQ(item->corner0.x, std::round(bb.x0));
                EXPECT_DOUBLE_EQ(item->corner1.x - item->corner0.x, im.cols);
                EXPECT_DOUBLE_EQ(item->corner1.y - item->corner0.y, im.rows);
                found = true;
            }
        }
        EXPECT_TRUE(found || im.rows == 0);
    }
}

INSTANTIATE_TEST_SUITE_P(Figures, StructureCase, ::testing::ValuesIn(structureCases()),
                         [](const ::testing::TestParamInfo<std::string>& i) { return i.param; });

}  // namespace
