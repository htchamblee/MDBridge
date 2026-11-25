#pragma once
// A drawn MDBridge figure (render/panels build_figure, figure_dump.hpp) compared with MDATK's figure
// record (tests/parity/make_figure_cases.py, mpl_walker.py) by its contents (owner decision D10;
// docs/specs/m4-plan.md section 3):
//
//   numbers (limits, line data, offsets, sizes, extents, norms)   1e-9 relative (1e-9 absolute
//                                                                 below 1); NaN equals NaN
//   colours                                                       within 1/255 per channel
//   strings, counts, flags, kinds                                 exact
//   layout boxes (axes, tight boxes, texts, titles, axis labels)  1.5 % of the figure size per edge
//   ticks and tick labels                                         exact whenever MDBridge's automatic
//                                                                 bin count equals Matplotlib's
//   images (the .npz: data and the RGBA shown) and colour-bar     data 1e-9, RGBA within 1 per
//   lookup tables                                                 channel, tables within 1/255
//
// What MDATK's colour bars hold inside (Matplotlib's mesh and divider artists) is the framework's,
// not a panel's: colour-bar axes are compared by their record (limits, ticks, label, norm, map,
// outline), not by their artists. The later positions compare the animated artists' changed fields.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "figure_records.hpp"
#include "mdbridge/render/figure_dump.hpp"

namespace mdbridge::test {

struct FigureComparison {
    double width = 1, height = 1;
    double worstBox = 0.0;
    std::string worstWhat;
    int ticksCompared = 0, ticksSkipped = 0;
    int numbers = 0;  // values compared

    static const core::Value* opt(const core::Value& v, std::string_view key) { return v.find(key); }
    static const core::Value& at(const core::Value& v, std::string_view key) {
        static const core::Value none;
        const core::Value* p = v.find(key);
        return p ? *p : none;
    }
    static double num(const core::Value& v) {
        if (const double* d = v.get<double>()) return *d;
        if (const std::int64_t* i = v.get<std::int64_t>()) return static_cast<double>(*i);
        if (const bool* b = v.get<bool>()) return *b ? 1.0 : 0.0;
        if (const core::Table* t = v.get<core::Table>(); t && t->size() == 1 && (*t)[0].first == "$float") {
            const std::string s = (*t)[0].second.stringOr();
            return s == "nan" ? std::nan("") : s == "inf" ? INFINITY : -INFINITY;
        }
        return std::nan("");
    }
    static std::string str(const core::Value& v) {
        if (const std::string* s = v.get<std::string>()) return *s;
        if (v.isNull()) return "None";
        if (const bool* b = v.get<bool>()) return *b ? "True" : "False";
        return std::to_string(num(v));
    }
    static bool isNumberList(const core::Value& v) {
        const core::Array* a = v.get<core::Array>();
        if (!a) return false;
        for (const core::Value& x : *a) {
            if (x.isArray() || x.isString()) return false;
        }
        return true;
    }

    bool close(double a, double b, double rel = 1e-9) const {
        if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
        if (std::isinf(a) || std::isinf(b)) return a == b;
        return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b));
    }

    // Recursive comparison of two record values: numbers within 1e-9, colours (4 numbers in [0, 1]
    // under a colour key) within 1/255, the rest exact.
    void value(const core::Value& mine, const core::Value& theirs, const std::string& what, bool colour = false) {
        if (theirs.isArray()) {
            if (!mine.isArray()) {
                ADD_FAILURE() << what << ": MDBridge has " << str(mine) << ", Matplotlib a list";
                return;
            }
            const core::Array& m = *mine.get<core::Array>();
            const core::Array& t = *theirs.get<core::Array>();
            if (m.size() != t.size()) {
                ADD_FAILURE() << what << ": " << m.size() << " values, Matplotlib " << t.size();
                return;
            }
            const bool isColour = colour && t.size() == 4 && isNumberList(theirs);
            for (std::size_t i = 0; i < t.size(); ++i) {
                if (isColour) {
                    ++numbers;
                    EXPECT_LE(std::fabs(num(m[i]) - num(t[i])), 1.0 / 255 + 1e-12) << what << "[" << i << "]";
                } else {
                    value(m[i], t[i], what + "[" + std::to_string(i) + "]", colour);
                }
            }
            return;
        }
        if (theirs.isTable() && !(theirs.get<core::Table>()->size() == 1 && (*theirs.get<core::Table>())[0].first == "$float")) {
            for (const auto& [k, v] : *theirs.get<core::Table>()) {
                const core::Value* m = mine.find(k);
                if (!m) {
                    ADD_FAILURE() << what << "." << k << " missing in MDBridge's record";
                    continue;
                }
                value(*m, v, what + "." + k, colour || k.find("color") != std::string::npos || k == "mfc" || k == "mec");
            }
            return;
        }
        const bool numeric = theirs.get<double>() || theirs.get<std::int64_t>() ||
                             (theirs.isTable() && !theirs.get<std::string>());
        if (numeric && !theirs.get<bool>()) {
            ++numbers;
            EXPECT_TRUE(close(num(mine), num(theirs))) << what << ": MDBridge " << num(mine) << ", Matplotlib " << num(theirs);
            return;
        }
        EXPECT_EQ(str(mine), str(theirs)) << what;
    }

    void box(const core::Value& mine, const core::Value& theirs, const std::string& what) {
        if (theirs.isNull()) return;
        if (mine.isNull()) {
            ADD_FAILURE() << what << ": no box in MDBridge's record";
            return;
        }
        const std::vector<double> m = numbersOf(mine), t = numbersOf(theirs);
        const double d = std::max({std::fabs(m[0] - t[0]) / width, std::fabs(m[2] - t[2]) / width,
                                   std::fabs(m[1] - t[1]) / height, std::fabs(m[3] - t[3]) / height});
        if (d > worstBox) {
            worstBox = d;
            worstWhat = what;
        }
        EXPECT_LE(d, 0.015) << what << ": MDBridge [" << m[0] << ", " << m[1] << ", " << m[2] << ", " << m[3]
                            << "] Matplotlib [" << t[0] << ", " << t[1] << ", " << t[2] << ", " << t[3] << "]";
    }

    // the record's fields `keys` of mine and theirs
    void fields(const core::Value& mine, const core::Value& theirs, std::initializer_list<const char*> keys,
                const std::string& what) {
        for (const char* k : keys) {
            const core::Value* t = theirs.find(k);
            if (!t) continue;
            const core::Value* m = mine.find(k);
            if (!m) {
                ADD_FAILURE() << what << "." << k << " missing in MDBridge's record";
                continue;
            }
            const std::string key(k);
            value(*m, *t, what + "." + key, key.find("color") != std::string::npos || key == "mfc" || key == "mec");
        }
    }

    void axis(const core::Value& mine, const core::Value& theirs, const std::string& what) {
        fields(at(mine, "label"), at(theirs, "label"), {"text", "color", "size"}, what + ".label");
        box(at(at(mine, "label"), "extent_px"), at(at(theirs, "label"), "extent_px"), what + " label");
        fields(mine, theirs, {"visible", "inverted"}, what);
        if (!at(theirs, "tick_style").isNull()) {
            fields(at(mine, "tick_style"), at(theirs, "tick_style"),
                   {"tick1", "tick2", "label1", "label2", "color", "length", "width", "labelcolor", "labelsize",
                    "rotation", "gridOn"},
                   what + ".tick_style");
        }
        const core::Value& tl = at(theirs, "locator");
        if (const core::Value* an = opt(tl, "auto_nbins")) {
            const core::Value* mn = opt(at(mine, "locator"), "auto_nbins");
            if (!mn || num(*mn) != num(*an)) {
                ++ticksSkipped;
                return;
            }
        }
        ++ticksCompared;
        // [loc, label, label1 shown, label2 shown]
        value(at(mine, "ticks"), at(theirs, "ticks"), what + ".ticks");
        EXPECT_EQ(str(at(mine, "offset_text")), str(at(theirs, "offset_text"))) << what;
        const std::string cls = str(at(tl, "class"));
        const std::string mcls = str(at(at(mine, "locator"), "class"));
        EXPECT_EQ(mcls == "FixedLocator", cls == "FixedLocator") << what << " locator " << mcls << " / " << cls;
        if (cls == "MaxNLocator") fields(at(mine, "locator"), tl, {"nbins", "integer", "min_n_ticks"}, what + ".locator");
    }

    void artist(const core::Value& mine, const core::Value& theirs, const std::string& what) {
        const std::string kind = str(at(theirs, "kind"));
        EXPECT_EQ(str(at(mine, "kind")), kind) << what;
        if (str(at(mine, "kind")) != kind) return;
        fields(mine, theirs, {"zorder", "visible", "alpha", "label"}, what);
        if (kind == "line") {
            fields(mine, theirs, {"x", "y", "color", "lw", "ls", "dashes", "marker", "ms", "mfc", "mec", "mew"}, what);
            // the walker tells transforms apart by equality: an axvline (x data, y axes) on an axes
            // whose x limits are 0..1 equals transAxes and is named "axes"
            const std::string mt = str(at(mine, "transform")), tt = str(at(theirs, "transform"));
            if (!(tt == "axes" && (mt == "xaxis" || mt == "yaxis"))) {
                EXPECT_EQ(mt, tt) << what << ".transform";
            }
        } else if (kind == "scatter") {
            fields(mine, theirs, {"offsets", "sizes", "facecolors", "edgecolors", "linewidths", "marker"}, what);
        } else if (kind == "segments") {
            fields(mine, theirs, {"segments", "colors", "linewidths", "capstyle"}, what);
        } else if (kind == "image") {
            fields(mine, theirs, {"shape", "extent", "origin", "cmap", "vmin", "vmax", "x_edges"}, what);
            EXPECT_EQ(opt(mine, "x_edges") != nullptr, opt(theirs, "x_edges") != nullptr) << what << " x_edges";
        } else if (kind == "text") {
            // (the walker names the coordinates of an axes text "data" when the axes' limits are 0..1)
            fields(mine, theirs, {"text", "ha", "va", "color", "size", "weight", "style", "rotation", "wrap"}, what);
            box(at(mine, "extent_px"), at(theirs, "extent_px"), what + " '" + str(at(theirs, "text")) + "'");
        } else if (kind == "legend") {
            fields(mine, theirs, {"loc"}, what);
            fields(at(mine, "frame"), at(theirs, "frame"), {"on", "facecolor", "edgecolor", "lw"}, what + ".frame");
            const core::Array& me = *at(mine, "entries").get<core::Array>();
            const core::Array& te = *at(theirs, "entries").get<core::Array>();
            ASSERT_EQ(me.size(), te.size()) << what << " entries";
            for (std::size_t i = 0; i < te.size(); ++i) {
                const std::string w = what + ".entries[" + std::to_string(i) + "]";
                fields(me[i], te[i], {"text", "color", "size"}, w);
                fields(at(me[i], "handle"), at(te[i], "handle"),
                       {"color", "lw", "ls", "dashes", "marker", "ms", "mfc", "mec", "mew"}, w + ".handle");
            }
            if (static_cast<int>(num(at(theirs, "loc"))) != 0) box(at(mine, "extent_px"), at(theirs, "extent_px"), what);
        }
    }

    // the whole first-position figure
    void figure(const core::Value& mine, const core::Value& theirs) {
        fields(mine, theirs, {"size_px", "dpi", "facecolor"}, "figure");
        if (!at(theirs, "suptitle").isNull()) {
            fields(at(mine, "suptitle"), at(theirs, "suptitle"), {"text", "color", "size", "ha", "va"}, "suptitle");
            box(at(at(mine, "suptitle"), "extent_px"), at(at(theirs, "suptitle"), "extent_px"), "suptitle");
        } else {
            EXPECT_TRUE(at(mine, "suptitle").isNull());
        }
        const core::Array& mt = *at(mine, "texts").get<core::Array>();
        const core::Array& tt = *at(theirs, "texts").get<core::Array>();
        ASSERT_EQ(mt.size(), tt.size());
        for (std::size_t i = 0; i < tt.size(); ++i) artist(mt[i], tt[i], "fig/text" + std::to_string(i));
        const core::Array& ma = *at(mine, "axes").get<core::Array>();
        const core::Array& ta = *at(theirs, "axes").get<core::Array>();
        ASSERT_EQ(ma.size(), ta.size()) << "axes";
        for (std::size_t k = 0; k < ta.size(); ++k) {
            const core::Value& m = ma[k];
            const core::Value& t = ta[k];
            const std::string id = "ax" + std::to_string(k) + " (" + str(at(t, "role")) + ")";
            fields(m, t, {"panel", "role", "facecolor", "visible", "axison", "xlim", "ylim", "zorder"}, id);
            box(at(m, "position_px"), at(t, "position_px"), id + " position");
            box(at(m, "tight_px"), at(t, "tight_px"), id + " tight box");
            for (const auto& [side, sp] : *at(t, "spines").get<core::Table>()) {
                if (side == "outline") continue;  // a colour bar's outline: colorbar.outline below
                fields(at(at(m, "spines"), side), sp, {"visible", "color", "lw"}, id + ".spines." + side);
            }
            for (const auto& [loc, title] : *at(t, "titles").get<core::Table>()) {
                const core::Value* mtitle = opt(at(m, "titles"), loc);
                ASSERT_NE(mtitle, nullptr) << id << " title " << loc;
                fields(*mtitle, title, {"text", "color", "size"}, id + ".title_" + loc);
                box(at(*mtitle, "extent_px"), at(title, "extent_px"), id + " title " + loc);
            }
            axis(at(m, "xaxis"), at(t, "xaxis"), id + ".xaxis");
            axis(at(m, "yaxis"), at(t, "yaxis"), id + ".yaxis");
            if (const core::Value* cb = opt(t, "colorbar")) {
                const core::Value* mcb = opt(m, "colorbar");
                ASSERT_NE(mcb, nullptr) << id;
                fields(*mcb, *cb, {"orientation", "vmin", "vmax", "label"}, id + ".colorbar");
                fields(at(*mcb, "cmap"), at(*cb, "cmap"), {"name", "N"}, id + ".colorbar.cmap");
                fields(at(*mcb, "outline"), at(*cb, "outline"), {"color", "lw"}, id + ".colorbar.outline");
                continue;  // the colour bar's own artists are Matplotlib's internals
            }
            std::map<std::string, const core::Value*> byId;
            for (const core::Value& a : *at(m, "artists").get<core::Array>()) byId[str(at(a, "id"))] = &a;
            std::set<std::string> seen;
            for (const core::Value& a : *at(t, "artists").get<core::Array>()) {
                const std::string aid = str(at(a, "id"));
                const std::string kind = str(at(a, "kind"));
                if (kind == "other") continue;
                const auto it = byId.find(aid);
                if (it == byId.end()) {
                    ADD_FAILURE() << aid << " (" << kind << ") missing in MDBridge's figure";
                    continue;
                }
                seen.insert(aid);
                artist(*it->second, a, aid);
            }
            for (const auto& [aid, a] : byId) EXPECT_TRUE(seen.count(aid)) << aid << " not in Matplotlib's figure";
        }
    }
};

// The walker's .npz archives (stored members, render_common.write_npz): name -> (shape, values as
// doubles) for the dtypes it writes ('<f8', '|u1', '<i8', '|b1')
struct NpyArray {
    std::vector<std::size_t> shape;
    std::vector<double> values;
};
inline std::map<std::string, NpyArray> readStoredNpz(const std::filesystem::path& path) {
    std::map<std::string, NpyArray> out;
    if (!std::filesystem::exists(path)) return out;
    const std::string zip = readFile(path);
    auto u16 = [&](std::size_t at) { return static_cast<std::size_t>(static_cast<unsigned char>(zip[at])) |
                                            static_cast<std::size_t>(static_cast<unsigned char>(zip[at + 1])) << 8; };
    auto u32 = [&](std::size_t at) { return u16(at) | u16(at + 2) << 16; };
    std::size_t at = 0;
    while (at + 30 <= zip.size() && u32(at) == 0x04034b50) {
        if (u16(at + 8) != 0) throw std::runtime_error(path.string() + ": a deflated member");
        const std::size_t size = u32(at + 18), nameLen = u16(at + 26), extraLen = u16(at + 28);
        std::string name = zip.substr(at + 30, nameLen);
        const std::size_t data = at + 30 + nameLen + extraLen;
        const std::string npy = zip.substr(data, size);
        at = data + size;
        const std::size_t headerLen = static_cast<std::size_t>(static_cast<unsigned char>(npy[8])) |
                                      static_cast<std::size_t>(static_cast<unsigned char>(npy[9])) << 8;
        const std::string header = npy.substr(10, headerLen);
        const std::string body = npy.substr(10 + headerLen);
        NpyArray a;
        const std::size_t s0 = header.find("'shape': (") + 10, s1 = header.find(')', s0);
        std::stringstream shape(header.substr(s0, s1 - s0));
        for (std::string dim; std::getline(shape, dim, ',');) {
            if (dim.find_first_not_of(' ') != std::string::npos) a.shape.push_back(std::stoul(dim));
        }
        const std::size_t d0 = header.find("'descr': '") + 10;
        const std::string descr = header.substr(d0, header.find('\'', d0) - d0);
        if (descr == "<f8") {
            a.values.resize(body.size() / 8);
            std::memcpy(a.values.data(), body.data(), a.values.size() * 8);
        } else if (descr == "<i8") {
            for (std::size_t i = 0; i + 8 <= body.size(); i += 8) {
                std::int64_t v;
                std::memcpy(&v, body.data() + i, 8);
                a.values.push_back(static_cast<double>(v));
            }
        } else if (descr == "|u1" || descr == "|b1") {
            for (char ch : body) a.values.push_back(static_cast<unsigned char>(ch));
        } else {
            throw std::runtime_error(path.string() + ": " + name + ": dtype " + descr);
        }
        if (name.size() > 4 && name.compare(name.size() - 4, 4, ".npy") == 0) name.resize(name.size() - 4);
        out[name] = std::move(a);
    }
    return out;
}

// The .npz arrays of position `pos` that the walker wrote for images and colour bars, compared with
// the figure's own: "p<pos>/ax<k>/<n>.data" (float64, NaN masked), ".rgba" (uint8), "p0/ax<k>.lut".
inline int compareArrays(const render::Figure& fig, const std::map<std::string, NpyArray>& npz, std::size_t pos) {
    int compared = 0;
    const std::string prefix = "p" + std::to_string(pos) + "/";
    for (const auto& [name, array] : npz) {
        if (name.rfind(prefix, 0) != 0) continue;
        const std::string rest = name.substr(prefix.size());
        const std::size_t slash = rest.find('/'), dot = rest.rfind('.');
        const std::size_t k = std::stoul(rest.substr(2, (slash == std::string::npos ? dot : slash) - 2));
        const std::string what = rest.substr(dot + 1);
        if (k >= fig.axes.size()) {
            ADD_FAILURE() << name << ": no axes " << k;
            continue;
        }
        const render::Axes& ax = fig.axes[k];
        if (what == "lut") {
            const render::Colormap& cmap = *ax.colorbar->cmap;
            const std::vector<double>& lut = array.values;
            if (lut.size() != static_cast<std::size_t>(cmap.size()) * 4) {
                ADD_FAILURE() << name << ": " << lut.size() / 4 << " entries, MDBridge " << cmap.size();
                continue;
            }
            for (int i = 0; i < cmap.size(); ++i) {
                const render::Rgba c = cmap.lut()[static_cast<std::size_t>(i)];
                const double mine[4] = {c.r, c.g, c.b, c.a};
                for (int ch = 0; ch < 4; ++ch) {
                    EXPECT_LE(std::fabs(mine[ch] - lut[static_cast<std::size_t>(4 * i + ch)]), 1.0 / 255 + 1e-12)
                        << name << " entry " << i;
                }
            }
            ++compared;
            continue;
        }
        const std::size_t n = std::stoul(rest.substr(slash + 1, dot - slash - 1));
        const render::Image* im = n < ax.artists.size() ? std::get_if<render::Image>(&ax.artists[n]) : nullptr;
        if (!im) {
            ADD_FAILURE() << name << ": not an image in MDBridge's figure";
            continue;
        }
        if (what == "data") {
            const std::vector<double>& data = array.values;
            EXPECT_EQ(array.shape, (std::vector<std::size_t>{static_cast<std::size_t>(im->rows), static_cast<std::size_t>(im->cols)}))
                << name;
            if (data.size() != im->values.size()) {
                ADD_FAILURE() << name << ": size " << im->values.size() << ", Matplotlib " << data.size();
                continue;
            }
            int bad = 0;
            for (std::size_t i = 0; i < data.size(); ++i) {
                const double a = im->values[i], b = data[i];
                const bool same = (std::isnan(a) && std::isnan(b)) || std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b));
                if (!same && bad++ < 5) ADD_FAILURE() << name << "[" << i << "]: MDBridge " << a << ", Matplotlib " << b;
            }
            EXPECT_EQ(bad, 0) << name;
        } else if (what == "rgba") {
            const std::vector<std::uint8_t> px = im->pixels();
            const std::vector<double>& want = array.values;
            if (px.size() != want.size()) {
                ADD_FAILURE() << name << ": size " << px.size() << ", Matplotlib " << want.size();
                continue;
            }
            int bad = 0;
            for (std::size_t i = 0; i < px.size(); ++i) {
                if (std::fabs(px[i] - want[i]) > 1.0 && bad++ < 5) {
                    ADD_FAILURE() << name << " byte " << i << ": MDBridge " << int(px[i]) << ", Matplotlib " << want[i];
                }
            }
            EXPECT_EQ(bad, 0) << name;
        }
        ++compared;
    }
    return compared;
}

}  // namespace mdbridge::test
