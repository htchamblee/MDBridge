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
#include <functional>
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
#include "mdbridge/render/example.hpp"
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
    Placeholder(const core::PanelSpec& spec, Context& context) : Panel(spec.type, spec, context) {}

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
                // row 0 at the top: the top edge at the box' rounded bottom plus the image height
                EXPECT_DOUBLE_EQ(item->corner0.y, b.fig->heightPx() - (std::round(bb.y0) + im.rows));
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

// ---- counterparts of MDATK's tests and the improvements -------------------------------------------

namespace {

core::Project structureProject() {
    core::Project p;
    p.layout.preset = "single";
    p.layout.panels = {core::PanelSpec{"structure", {}, "", {}}};
    return p;
}

core::Value record(const std::string& name) { return test::untag(test::readJson(test::kFigureCases / (name + ".json"))); }

// a case's project and input, built at the case's size (the framework's own layout)
Built buildRecord(const std::string& name, const std::function<void(core::Project&)>& edit = {}) {
    const Value rec = record(name);
    core::Project p = test::caseProject(rec);
    if (edit) edit(p);
    const std::vector<double> size = test::numbersOf(at(rec, "size"));
    return buildFigure(test::loadRenderInput(str(at(rec, "input"))), p, static_cast<int>(size[0]),
                       static_cast<int>(size[1]), num(at(rec, "dpi")));
}

bool contains(const std::vector<ArtistRef>& refs, const ArtistRef& r) {
    return std::find(refs.begin(), refs.end(), r) != refs.end();
}

std::size_t indexOf(const std::vector<ArtistRef>& refs, const ArtistRef& r) {
    return static_cast<std::size_t>(std::find(refs.begin(), refs.end(), r) - refs.begin());
}

}  // namespace

TEST(StructurePanel, ZoomAndCentre) {
    // test_seqpanels.py test_structure_zoom_and_centre
    core::Project p = structureProject();
    const auto make = [&] { return buildFigure(exampleDataset(p), p, 800, 800, 100.0); };
    Built a = make();
    const double w1 = a.fig->axes[a.structure->mainAxes()].xaxis.limits.hi - a.fig->axes[a.structure->mainAxes()].xaxis.limits.lo;
    p.visual.zoom = 2.0;
    Built b = make();
    const Range xl = b.fig->axes[b.structure->mainAxes()].xaxis.limits;
    EXPECT_LT(xl.hi - xl.lo, 0.6 * w1);
    p.visual.center_on = "origin";
    Built c = make();
    for (double v : c.structure->center()) EXPECT_EQ(v, 0.0);
    p.visual.center_on = "A:6";
    Built d = make();
    const std::vector<double> coords = d.dataset->result("trace")->array("coords")->toDouble();
    for (std::size_t k = 0; k < 3; ++k) EXPECT_DOUBLE_EQ(d.structure->center()[k], coords[5 * 3 + k]);
    p.visual.center_on = "Z:99";  // not found: the model centre, with a note
    resetWarnOnce();
    LogCapture log;
    Built e = make();
    EXPECT_TRUE(e.structure->ok());
    EXPECT_NE(e.structure->viewNote().find("centre"), std::string::npos);
    ASSERT_EQ(log.lines.size(), 1u);
    EXPECT_EQ(log.lines[0], "WARNING Structure view \xe2\x80\xba Orientation \xe2\x80\xba Centre on: no residue of the "
                            "structure view matches 'Z:99'; using the model centre");
}

TEST(StructurePanel, FollowedCentreMovesTheModelWithTheResidue) {
    // test_seqpanels.py test_followed_centre_moves_the_model_with_the_residue
    core::Project p = structureProject();
    core::Dataset ds = exampleDataset(p);
    core::AnalysisResult tr = *ds.result("trace");
    const core::NumArray& c = *tr.array("coords");
    std::vector<double> coords = c.toDouble();
    const std::size_t frames = c.shape()[0], nRes = c.shape()[1];
    for (std::size_t f = 0; f < frames; ++f) {
        for (std::size_t r = 0; r < nRes; ++r) coords[(f * nRes + r) * 3] += 0.5 * static_cast<double>(f);  // drift along x
    }
    tr.setArray("coords", core::NumArray(coords, c.shape()));
    ds.add(tr);
    p.visual.center_on = "A:3";
    p.visual.center_follow = true;
    Built b = buildFigure(std::move(ds), p, 600, 600, 100.0);
    StructurePanel& st = *b.structure;
    for (std::size_t f = 0; f < frames; ++f) {
        for (std::size_t k = 0; k < 3; ++k) EXPECT_NEAR(st.view(f)[2][k], st.view(0)[2][k], 1e-9);  // stays in the middle
    }
    EXPECT_NEAR(st.frameShift().back()[0], 0.5 * static_cast<double>(frames - 1), 1e-9);
}

TEST(StructurePanel, ShadowFollowsTheLayoutAndAFollowedCentre) {
    // test_density.py: the shadow window follows the axes limits after layout (a wide panel), and
    // the shadow moves with a followed centre residue (animated, not a rotation artist)
    Built b = buildRecord("density_follow");
    StructurePanel& st = *b.structure;
    ASSERT_TRUE(st.hasDensity());
    const Axes& ax = b.fig->axes[st.mainAxes()];
    const Range xl = ax.xaxis.limits, yl = ax.yaxis.limits;
    EXPECT_GT(xl.hi - xl.lo, 1.5 * (yl.hi - yl.lo));
    const std::size_t di = *st.artists().density;
    const ArtistRef density{st.mainAxes(), di, false};
    EXPECT_TRUE(contains(st.animatedArtists(), density));
    EXPECT_TRUE(st.rotationArtists().empty());
    const std::array<double, 4> base = std::get<Image>(ax.artists[di]).extent;
    const double pad = 0.15 * (xl.hi - xl.lo) + 3.0;
    EXPECT_LE(xl.lo - pad, base[0]);
    EXPECT_LE(base[1], xl.hi + pad);
    // the frame where the residue moved most
    std::size_t pos = 0;
    double best = 0.0;
    for (std::size_t f = 0; f < st.frameShift().size(); ++f) {
        const auto& s = st.frameShift()[f];
        const double n = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        if (n > best) best = n, pos = f;
    }
    ASSERT_GT(best, 0.05);
    b.update(pos);
    const Mat3& v = st.viewMatrix();
    const auto& s = st.frameShift()[pos];
    const double dx = s[0] * v[0][0] + s[1] * v[0][1] + s[2] * v[0][2], dy = s[0] * v[1][0] + s[1] * v[1][1] + s[2] * v[1][2];
    const std::array<double, 4> moved = std::get<Image>(b.fig->axes[st.mainAxes()].artists[di]).extent;
    EXPECT_NEAR(moved[0], base[0] - dx, 1e-6);
    EXPECT_NEAR(moved[1], base[1] - dx, 1e-6);
    EXPECT_NEAR(moved[2], base[2] - dy, 1e-6);
    EXPECT_NEAR(moved[3], base[3] - dy, 1e-6);
}

TEST(StructurePanel, DensityPlacementsDegrade) {
    // test_density.py: a map that cannot be read, an unusable placement or none: a note, the panel
    // stays (no shadow)
    resetWarnOnce();
    {
        LogCapture log;
        Built b = buildRecord("density_behind", [](core::Project& p) { p.density.map_path = test::kData.string() + "/nowhere.mrc"; });
        EXPECT_TRUE(b.structure->ok());
        EXPECT_FALSE(b.structure->hasDensity());
        EXPECT_EQ(b.structure->densityNote().rfind("density map: ", 0), 0u) << b.structure->densityNote();
        ASSERT_FALSE(log.lines.empty());
        EXPECT_EQ(log.lines[0].rfind("WARNING density map not drawn: ", 0), 0u);
    }
    {
        Built b = buildRecord("density_behind", [](core::Project& p) {
            p.density.transform = {1, 0, 0, 0, 1, 0, 0, 0, 1, std::numeric_limits<double>::infinity(), 0, 0};
        });
        EXPECT_TRUE(b.structure->ok());
        EXPECT_EQ(b.structure->densityNote(), "density map: the transform contains non-finite numbers");
    }
    {
        Built b = buildRecord("density_behind", [](core::Project& p) { p.density.transform.clear(); });
        EXPECT_EQ(b.structure->densityNote(), "density map: not fitted yet");
        EXPECT_FALSE(b.structure->artists().density.has_value());
    }
}

TEST(StructurePanel, DensityLevelFromFiniteVoxels) {
    // RN-33: a map with NaN voxels still shows its shadow (MDATK: the level became NaN)
    const Built clean = buildRecord("density_behind");
    io::DensityMap map = io::readMap(test::kData / "render/inputs/demo_dry_map.mrc");
    const double level = io::mapLevel(map, 1.0);
    EXPECT_EQ(finiteMapLevel(map, 1.0), level);  // all finite: MDATK's value
    map.data[0] = std::nanf("");
    map.data[1] = std::numeric_limits<float>::infinity();
    EXPECT_TRUE(std::isnan(io::mapLevel(map, 1.0)));
    const double finite = finiteMapLevel(map, 1.0);
    EXPECT_TRUE(std::isfinite(finite));
    EXPECT_NEAR(finite, level, 1e-3 * std::fabs(level) + 1e-6);  // two corner voxels of a large map
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "mdbridge_w2b_nan_map.mrc";
    io::writeMap(path, map);
    Built b = buildRecord("density_behind", [&](core::Project& p) { p.density.map_path = path.string(); });
    ASSERT_TRUE(b.structure->hasDensity());
    const Image& mine = std::get<Image>(b.fig->axes[b.structure->mainAxes()].artists[*b.structure->artists().density]);
    const Image& ref = std::get<Image>(clean.fig->axes[clean.structure->mainAxes()].artists[*clean.structure->artists().density]);
    EXPECT_TRUE(mine.visible);
    EXPECT_EQ(mine.rows, ref.rows);
    EXPECT_EQ(mine.cols, ref.cols);
    int alpha = 0;
    for (std::size_t i = 3; i < mine.rgba.size(); i += 4) alpha += mine.rgba[i] > 0 ? 1 : 0;
    EXPECT_GT(alpha, 0);
    std::filesystem::remove(path);
    EXPECT_TRUE(std::isnan(finiteMapLevel(io::DensityMap{{1, 1, 1}, {std::nanf("")}, {}, {1, 1, 1}, {}}, 1.0)));
}

TEST(StructurePanel, RotationNudgeAndFastQuality) {
    // set_rotation re-projects the model and the shadow; nudge_density moves the map on screen;
    // fast quality draws a coarser shadow and a half-resolution 3D layer at the full size
    Built b = buildRecord("density_front_rep");
    StructurePanel& st = *b.structure;
    Figure& fig = *b.fig;
    const std::size_t di = *st.artists().density;
    const auto extent = [&] { return std::get<Image>(fig.axes[st.mainAxes()].artists[di]).extent; };
    const std::array<double, 4> e0 = extent();
    const Vec3d v0 = st.view(0)[0];
    const std::vector<double> angles{0.0, 90.0, 0.0};
    st.setRotation(fig, angles);
    EXPECT_NE(st.view(0)[0], v0);
    EXPECT_NE(extent(), e0);
    st.setRotation(fig, b.project->visual.rotation);
    const std::array<double, 4> e1 = extent();
    for (std::size_t k = 0; k < 4; ++k) EXPECT_NEAR(e1[k], e0[k], 1e-9);
    const core::RigidTransform before = *st.densityPlacement();
    const std::optional<std::vector<double>> placed = st.nudgeDensity(fig, {3.0, -2.0});
    ASSERT_TRUE(placed.has_value());
    EXPECT_EQ(placed->size(), 12u);
    const std::array<double, 4> e2 = extent();
    // the shadow moves with the map: by the shift on screen (to the voxel step of its grid)
    EXPECT_NEAR(0.5 * (e2[0] + e2[1]) - 0.5 * (e1[0] + e1[1]), 3.0, 2.0);
    EXPECT_NEAR(0.5 * (e2[2] + e2[3]) - 0.5 * (e1[2] + e1[3]), -2.0, 2.0);
    // a rotation about the map's density centre keeps the centre where it is
    const Mat3 r = rotationZ(30.0);
    ASSERT_TRUE(st.nudgeDensity(fig, {0.0, 0.0}, r).has_value());
    EXPECT_NE(st.densityPlacement()->rotation, before.rotation);
    // fast quality
    const int cols = std::get<Image>(fig.axes[st.mainAxes()].artists[di]).cols;
    st.setFast(true);
    st.setRotation(fig, b.project->visual.rotation);
    b.update(0);
    EXPECT_EQ(st.stageQuality(), "fast");
    EXPECT_LT(std::get<Image>(fig.axes[st.mainAxes()].artists[di]).cols, cols + 1);
    const Image& stage = std::get<Image>(fig.axes[st.mainAxes()].artists[*st.artists().stage]);
    const Box bb = fig.axesBox(st.mainAxes());
    EXPECT_EQ(stage.cols, static_cast<int>(std::lround(bb.width())));
    EXPECT_EQ(stage.rows, static_cast<int>(std::lround(bb.height())));
}

TEST(StructurePanel, FrontShadowIsAnimatedAboveThe3dLayer) {
    // test_molecule_render.py test_front_density_shadow_is_animated_above_the_3d_layer
    Built b = buildRecord("density_front_rep");
    StructurePanel& st = *b.structure;
    const ArtistRef density{st.mainAxes(), *st.artists().density, false};
    const ArtistRef stage{st.mainAxes(), *st.artists().stage, false};
    std::vector<ArtistRef> arts = st.animatedArtists();
    EXPECT_GT(indexOf(arts, density), indexOf(arts, stage));
    EXPECT_TRUE(st.rotationArtists().empty());
    Built behind = buildRecord("density_front_rep", [](core::Project& p) { p.density.front = false; });
    StructurePanel& sb = *behind.structure;
    const ArtistRef d2{sb.mainAxes(), *sb.artists().density, false};
    EXPECT_FALSE(contains(sb.animatedArtists(), d2));
    EXPECT_EQ(sb.rotationArtists(), std::vector<ArtistRef>{d2});
    // the frame loop's layers: the 3D layer and the front shadow are drawn per frame
    markAnimated(*b.fig, st.animatedArtists());
    const DisplayList dyn = buildDisplayList(*b.fig, fonts(), DisplayLayer::Animated);
    int images = 0;
    for (const DisplayItem& d : dyn.items) images += std::holds_alternative<ImageItem>(d.item) ? 1 : 0;
    EXPECT_EQ(images, 2);
}

namespace {

class FailingRenderer final : public StageRenderer {
public:
    molecule::RgbaImage render(const std::vector<molecule::Colored>&, const molecule::Viewport&, const molecule::Shading&,
                               molecule::RenderQuality) const override {
        throw std::runtime_error("no luck");
    }
    std::string name() const override { return "FailingRenderer"; }
};

class FailingFactory final : public StageRendererFactory {
public:
    StageRendererChoice choose(std::string_view) const override { return {std::make_shared<FailingRenderer>(), {}}; }
};

}  // namespace

TEST(StructurePanel, FailedRenderFallsBackToThe2dTrace) {
    // test_molecule_render.py test_failed_3d_render_falls_back_to_the_2d_trace
    LogCapture log;
    setStageRendererFactory(std::make_shared<FailingFactory>());
    Built b = buildRecord("rep_cartoon");
    setStageRendererFactory(nullptr);
    StructurePanel& st = *b.structure;
    b.update(2);
    b.update(3);
    EXPECT_TRUE(st.stageFailed());
    const Axes& ax = b.fig->axes[st.mainAxes()];
    EXPECT_EQ(std::get<Scatter>(ax.artists[st.artists().points]).offsets.size(), st.view(0).size());
    EXPECT_FALSE(std::get<Scatter>(ax.artists[st.artists().ions]).offsets.empty());  // ions in 2D too
    EXPECT_FALSE(std::get<Segments>(ax.artists[st.artists().lines]).segments.empty());
    int warnings = 0;
    for (const std::string& l : log.lines) {
        warnings += l == "WARNING 3D representation not drawn (2D trace shown instead): no luck" ? 1 : 0;
    }
    EXPECT_EQ(warnings, 1);  // once, not per frame
    // and back to 3D when a render succeeds again
    Built ok = buildRecord("rep_cartoon");
    EXPECT_FALSE(ok.structure->stageFailed());
    EXPECT_TRUE(std::get<Scatter>(ok.fig->axes[ok.structure->mainAxes()].artists[ok.structure->artists().points]).offsets.empty());
}

TEST(StructurePanel, NonFiniteFrameIsLeftOut) {
    // RN-03: a non-finite trace point (frame 3 of demo_dry_nan): finite extent, one warning, the
    // point and its links not drawn in that frame
    resetWarnOnce();
    LogCapture log;
    Built b = buildRecord("structure_nan");
    StructurePanel& st = *b.structure;
    EXPECT_TRUE(std::isfinite(st.radius()));
    b.update(3);
    const Axes& ax = b.fig->axes[st.mainAxes()];
    const Scatter& pts = std::get<Scatter>(ax.artists[st.artists().points]);
    EXPECT_EQ(pts.offsets.size(), st.view(3).size() - 1);
    for (const Point& p : pts.offsets) EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y));
    for (const auto& seg : std::get<Segments>(ax.artists[st.artists().lines]).segments) {
        for (const Point& p : seg) EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y));
    }
    for (double s : pts.sizes) EXPECT_TRUE(std::isfinite(s));
    ASSERT_EQ(log.lines.size(), 1u);
    EXPECT_EQ(log.lines[0], "WARNING structure view: non-finite coordinates in frame(s) 3; left out of the view "
                            "(those residues are not drawn there)");
    b.update(4);
    EXPECT_EQ(std::get<Scatter>(b.fig->axes[st.mainAxes()].artists[st.artists().points]).offsets.size(), st.view(4).size());
}

TEST(StructurePanel, IonAndChainLegend) {
    // RN-08: the ions' element; RN-09: chains by name (D11)
    EXPECT_EQ(ionLegendLabel({"MG25", "MG26"}), "Mg\xc2\xb2\xe2\x81\xba");
    EXPECT_EQ(ionLegendLabel({"MG225"}), "Mg\xc2\xb2\xe2\x81\xba");  // resname MG2
    EXPECT_EQ(ionLegendLabel({"NA101", "Na+102", "SOD7"}), "Na\xe2\x81\xba");
    EXPECT_EQ(ionLegendLabel({"K1"}), "K\xe2\x81\xba");
    EXPECT_EQ(ionLegendLabel({"CA5"}), "Ca\xc2\xb2\xe2\x81\xba");
    EXPECT_EQ(ionLegendLabel({"MG1", "K2"}), "Ions");
    EXPECT_EQ(ionLegendLabel({"XYZ1"}), "Ions");
    EXPECT_EQ(ionLegendLabel({}), "Ions");
    Built b = buildRecord("structure_chain");
    const Legend& L = *b.fig->axes[b.structure->mainAxes()].legend;
    ASSERT_GE(L.entries.size(), 2u);
    EXPECT_EQ(L.entries[0].text, "Chain A");
    EXPECT_EQ(L.entries[1].text, "Chain B");
    EXPECT_EQ(L.entries.back().text, "Mg\xc2\xb2\xe2\x81\xba");
}

TEST(StructurePanel, GpuRendererNote) {
    // RN-16: "vtk" / "vulkan" mean the GPU renderer; without one the built-in renderer with a note
    EXPECT_EQ(defaultStageRenderer("vtk").note, kNoGpuRendererNote);
    EXPECT_EQ(defaultStageRenderer("vulkan").note, kNoGpuRendererNote);
    EXPECT_TRUE(defaultStageRenderer("auto").note.empty());
    EXPECT_EQ(defaultStageRenderer("software").renderer->name(), "SoftwareRenderer");
}

TEST(StructurePanel, MissingResultsShowNoData) {
    // MDATK's KeyError: "No data: trace" (and "No data: atoms" for representations without atoms)
    core::Project p = structureProject();
    core::Dataset empty({0, 1}, {0.0, 10.0});
    Built b = buildFigure(std::move(empty), p, 640, 360, 60.0);
    EXPECT_FALSE(b.structure->ok());
    const Text& t = std::get<Text>(b.fig->axes[b.structure->mainAxes()].artists.at(0));
    EXPECT_EQ(t.text, "No data: trace");
    EXPECT_TRUE(b.structure->animatedArtists().empty());
    Built c = buildRecord("rep_cartoon");
    core::Dataset noAtoms(c.dataset->frames(), c.dataset->timePs(), {}, c.dataset->timeUnit());
    for (const core::AnalysisResult& r : c.dataset->results()) {
        if (r.name != "atoms") noAtoms.add(r);
    }
    Built d = buildFigure(std::move(noAtoms), *c.project, 640, 360, 60.0);
    EXPECT_FALSE(d.structure->ok());
    EXPECT_EQ(std::get<Text>(d.fig->axes[d.structure->mainAxes()].artists.at(0)).text, "No data: atoms");
    // proximity colours without a proximity result; the frame stays (RN-42), built as build_figure does
    core::Project q = *c.project;
    q.visual.representations.clear();
    q.visual.color_mode = "proximity";
    core::Dataset noProx(c.dataset->frames(), c.dataset->timePs(), {}, c.dataset->timeUnit());
    noProx.add(*c.dataset->result("trace"));
    Context context(noProx, q);
    Figure fig(640, 360, 60.0);
    std::unique_ptr<Panel> panel = createPanel(fig, {0, 1, 0, 1}, q.layout.panels.at(0), context);
    EXPECT_FALSE(panel->ok());
    const Axes& ax = fig.axes[panel->mainAxes()];
    EXPECT_EQ(std::get<Text>(ax.artists.at(0)).text, "No data: proximity");
    EXPECT_TRUE(ax.spines[kLeft].visible);
    EXPECT_TRUE(ax.spines[kBottom].visible);
    EXPECT_FALSE(ax.spines[kTop].visible);
}

TEST(StructurePanel, RegisteredAndBuiltByBuildFigure) {
    // the registry knows the structure panel; build_figure of a structure case gives the same view
    const std::vector<std::string_view> types = registeredPanelTypes();
    EXPECT_NE(std::find(types.begin(), types.end(), "structure"), types.end());
    const Value rec = record("structure_tags");
    const core::Dataset ds = test::loadRenderInput(str(at(rec, "input")));
    const core::Project p = test::caseProject(rec);
    BuiltFigure built = test::buildCase(rec, ds, p);
    ASSERT_EQ(built.panels.size(), 1u);
    EXPECT_EQ(built.panels[0]->type(), "structure");
    EXPECT_TRUE(built.panels[0]->ok());
    built.update(12, fonts());
    EXPECT_FALSE(built.animatedArtists().empty());
}
