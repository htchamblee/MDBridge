// The glob-pattern translation of fnmatch.fnmatchcase (Lib/fnmatch.py)
// is ported from CPython (Python Software Foundation License Version 2; full text in
// packaging/licenses/Python-LICENSE.txt, see THIRD_PARTY_NOTICES.md):
// Copyright (c) 2001-2024 Python Software Foundation; All Rights Reserved

#include "mdbridge/core/atom_selection.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

#include "selection_tables.inc"

using Mask = std::vector<std::uint8_t>;
using Token = std::optional<std::string>;

std::string show(const Token& token) { return token ? *token : "None"; }

// ---- Python values ---------------------------------------------------------------------------

// A Python int compared with int64 values: `overflow` > 0 above every int64, < 0 below every one.
struct PyInt {
    int overflow = 0;
    std::int64_t value = 0;
};

bool equal(std::int64_t v, const PyInt& x) { return x.overflow == 0 && v == x.value; }
bool atLeast(std::int64_t v, const PyInt& x) {
    return x.overflow == 0 ? v >= x.value : x.overflow < 0;
}
bool atMost(std::int64_t v, const PyInt& x) {
    return x.overflow == 0 ? v <= x.value : x.overflow > 0;
}

// int(text); throws ValueError with Python's message
PyInt pyInt(const std::string& text) {
    try {
        return {0, pythonInt(text)};
    } catch (const std::overflow_error&) {
        const std::string stripped = pythonStrip(text);
        return {!stripped.empty() && stripped[0] == '-' ? -1 : 1, 0};
    }
}

PyInt pyInt(const Token& token) {
    if (!token) {
        throw AtomSelectionError(
            "int() argument must be a string, a bytes-like object or a real number, not "
            "'NoneType'");
    }
    return pyInt(*token);
}

double pyFloat(const Token& token) {
    if (!token) {
        throw AtomSelectionError(
            "float() argument must be a string or a real number, not 'NoneType'");
    }
    return pythonFloat(*token);
}

// numpy's cast of a Python float to float32 (overflow to infinity)
float toFloat32(double x) {
    constexpr double kLimit = 3.4028235677973366e38;  // halfway between FLT_MAX and 2^128
    if (std::isnan(x)) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    if (x >= kLimit) {
        return std::numeric_limits<float>::infinity();
    }
    if (x <= -kLimit) {
        return -std::numeric_limits<float>::infinity();
    }
    return static_cast<float>(x);
}

// np.isclose(x, y) with the default tolerances, y a Python float
bool isClose(double x, double y) {
    return (std::fabs(x - y) <= 1e-8 + 1e-5 * std::fabs(y) && std::isfinite(y)) || x == y;
}

// np.isclose(float32 array, Python float): the difference and the comparisons in float32
bool isClose32(float x, double y) {
    const float yf = toFloat32(y);
    const float threshold = toFloat32(1e-8 + 1e-5 * std::fabs(y));
    return (std::fabs(x - yf) <= threshold && std::isfinite(y)) || x == yf;
}

// np.isclose(float32 array, float32 scalar): everything in float32
bool isClose32(float x, float y) {
    const float threshold = 1e-8f + 1e-5f * std::fabs(y);
    return (std::fabs(x - y) <= threshold && std::isfinite(y)) || x == y;
}

// ---- regular expressions of the parser (re.match semantics: backtracking, prefix) -------------

struct Match {
    const std::u32string& s;
    std::array<std::pair<std::size_t, std::size_t>, 4> groups{};
};
using Next = std::function<bool(std::size_t)>;
using Rx = std::function<bool(Match&, std::size_t, const Next&)>;

Rx charIf(std::function<bool(char32_t)> predicate) {
    return [predicate](Match& m, std::size_t pos, const Next& next) {
        return pos < m.s.size() && predicate(m.s[pos]) && next(pos + 1);
    };
}

Rx literal(std::u32string text) {
    return [text](Match& m, std::size_t pos, const Next& next) {
        return m.s.compare(pos, text.size(), text) == 0 && next(pos + text.size());
    };
}

Rx optional(Rx rx) {
    return
        [rx](Match& m, std::size_t pos, const Next& next) { return rx(m, pos, next) || next(pos); };
}

Rx star(Rx rx) {
    auto self = std::make_shared<Rx>();
    *self = [rx, weak = std::weak_ptr<Rx>(self)](Match& m, std::size_t pos, const Next& next) {
        const auto again = weak.lock();
        return rx(m, pos, [&](std::size_t q) { return q != pos && (*again)(m, q, next); }) ||
               next(pos);
    };
    return [self](Match& m, std::size_t pos, const Next& next) { return (*self)(m, pos, next); };
}

Rx sequence(std::vector<Rx> parts) {
    return [parts](Match& m, std::size_t pos, const Next& next) {
        std::function<bool(std::size_t, std::size_t)> step = [&](std::size_t k, std::size_t at) {
            if (k == parts.size()) {
                return next(at);
            }
            return parts[k](m, at, [&](std::size_t q) { return step(k + 1, q); });
        };
        return step(0, pos);
    };
}

Rx plus(const Rx& rx) { return sequence({rx, star(rx)}); }

Rx either(Rx a, Rx b) {
    return [a, b](Match& m, std::size_t pos, const Next& next) {
        return a(m, pos, next) || b(m, pos, next);
    };
}

Rx group(std::size_t index, Rx rx) {
    return [index, rx](Match& m, std::size_t pos, const Next& next) {
        return rx(m, pos, [&](std::size_t q) {
            const auto saved = m.groups[index];
            m.groups[index] = {pos, q};
            if (next(q)) {
                return true;
            }
            m.groups[index] = saved;
            return false;
        });
    };
}

Rx end() {
    return
        [](Match& m, std::size_t pos, const Next& next) { return pos == m.s.size() && next(pos); };
}

Rx ch(char32_t c) {
    return charIf([c](char32_t x) { return x == c; });
}

const Rx& digit() {
    static const Rx rx = charIf(isPythonDigit);
    return rx;
}

const Rx& intPattern() {  // -?\d+
    static const Rx rx = sequence({optional(ch('-')), plus(digit())});
    return rx;
}

const Rx& floatPattern() {  // -?\d*\.?\d*(?:e[-+]?\d+)?
    static const Rx rx =
        sequence({optional(ch('-')), star(digit()), optional(ch('.')), star(digit()),
                  optional(sequence(
                      {ch('e'), optional(charIf([](char32_t c) { return c == '-' || c == '+'; })),
                       plus(digit())}))});
    return rx;
}

const Rx& rangeSeparator() {  // \s*(?:[:-]| to )\s*
    static const Rx rx =
        sequence({star(charIf(isPythonSpace)),
                  either(charIf([](char32_t c) { return c == ':' || c == '-'; }), literal(U" to ")),
                  star(charIf(isPythonSpace))});
    return rx;
}

const Rx& wordOptional() {  // \w?
    static const Rx rx = optional(charIf(isPythonWordCharacter));
    return rx;
}

// re.match(pattern, text): the groups as UTF-8, std::nullopt without a match
std::optional<std::vector<std::string>> reMatch(const Rx& rx, std::size_t groupCount,
                                                const std::string& text) {
    const std::u32string s = decodeUtf8(text);
    Match m{s};
    if (!rx(m, 0, [](std::size_t) { return true; })) {
        return std::nullopt;
    }
    std::vector<std::string> out;
    for (std::size_t i = 0; i < groupCount; ++i) {
        out.push_back(encodeUtf8(std::u32string_view(s).substr(
            m.groups[i].first, m.groups[i].second - m.groups[i].first)));
    }
    return out;
}

// ---- fnmatch.fnmatchcase ---------------------------------------------------------------------

struct GlobItem {
    enum class Kind { Star, Any, Literal, Set, Never } kind = Kind::Literal;
    char32_t ch = 0;
    bool negated = false;
    std::vector<std::pair<char32_t, char32_t>> ranges;

    bool matches(char32_t c) const {
        switch (kind) {
        case Kind::Any:
            return true;
        case Kind::Literal:
            return c == ch;
        case Kind::Never:
            return false;
        case Kind::Set: {
            bool in = false;
            for (const auto& [lo, hi] : ranges) {
                in = in || (c >= lo && c <= hi);
            }
            return in != negated;
        }
        case Kind::Star:
            break;
        }
        return false;
    }
};

std::u32string escapeSetText(const std::u32string& text, bool hyphens) {
    std::u32string out;
    for (const char32_t c : text) {
        if (c == U'\\' || (hyphens && c == U'-')) {
            out += U'\\';
        }
        out += c;
    }
    return out;
}

// The character class Python's re builds from fnmatch's "[...]" translation.
GlobItem setItem(const std::u32string& pat, std::size_t i, std::size_t j) {
    std::u32string stuff;
    const std::u32string inner = pat.substr(i, j - i);
    if (inner.find(U'-') == std::u32string::npos) {
        stuff = escapeSetText(inner, false);
    } else {
        std::vector<std::u32string> chunks;
        std::size_t k = pat[i] == U'!' ? i + 2 : i + 1;
        for (;;) {
            const std::size_t dash = k < j ? pat.find(U'-', k) : std::u32string::npos;
            if (dash == std::u32string::npos || dash >= j) {
                break;
            }
            chunks.push_back(pat.substr(i, dash - i));
            i = dash + 1;
            k = dash + 3;
        }
        const std::u32string chunk = pat.substr(i, j - i);
        if (!chunk.empty()) {
            chunks.push_back(chunk);
        } else {
            chunks.back() += U'-';
        }
        for (std::size_t m = chunks.size() - 1; m > 0; --m) {
            if (chunks[m - 1].empty() || chunks[m].empty()) {
                throw AtomSelectionError("string index out of range");
            }
            if (chunks[m - 1].back() > chunks[m].front()) {
                chunks[m - 1] =
                    chunks[m - 1].substr(0, chunks[m - 1].size() - 1) + chunks[m].substr(1);
                chunks.erase(chunks.begin() + static_cast<std::ptrdiff_t>(m));
            }
        }
        for (std::size_t m = 0; m < chunks.size(); ++m) {
            stuff += (m ? U"-" : U"") + escapeSetText(chunks[m], true);
        }
    }
    GlobItem item;
    if (stuff.empty()) {
        item.kind = GlobItem::Kind::Never;
        return item;
    }
    if (stuff == U"!") {
        item.kind = GlobItem::Kind::Any;
        return item;
    }
    item.kind = GlobItem::Kind::Set;
    std::size_t pos = 0;
    if (stuff[0] == U'!') {
        item.negated = true;
        pos = 1;
    }
    // the class body as sre reads it: atoms (\x escapes one character), a-b ranges
    auto atom = [&](std::size_t& at) {
        char32_t c = stuff[at++];
        if (c == U'\\' && at < stuff.size()) {
            c = stuff[at++];
        }
        return c;
    };
    while (pos < stuff.size()) {
        const char32_t lo = atom(pos);
        if (pos + 1 < stuff.size() && stuff[pos] == U'-') {
            ++pos;
            item.ranges.emplace_back(lo, atom(pos));
        } else {
            item.ranges.emplace_back(lo, lo);
        }
    }
    return item;
}

std::vector<GlobItem> compileGlob(const std::u32string& pat) {
    std::vector<GlobItem> out;
    std::size_t i = 0;
    const std::size_t n = pat.size();
    while (i < n) {
        const char32_t c = pat[i++];
        GlobItem item;
        if (c == U'*') {
            if (out.empty() || out.back().kind != GlobItem::Kind::Star) {
                item.kind = GlobItem::Kind::Star;
                out.push_back(item);
            }
            continue;
        }
        if (c == U'?') {
            item.kind = GlobItem::Kind::Any;
        } else if (c == U'[') {
            std::size_t j = i;
            if (j < n && pat[j] == U'!') {
                ++j;
            }
            if (j < n && pat[j] == U']') {
                ++j;
            }
            while (j < n && pat[j] != U']') {
                ++j;
            }
            if (j >= n) {
                item.ch = U'[';
            } else {
                item = setItem(pat, i, j);
                i = j + 1;
            }
        } else {
            item.ch = c;
        }
        out.push_back(item);
    }
    return out;
}

bool globMatch(const std::vector<GlobItem>& pattern, const std::u32string& text) {
    std::size_t p = 0, t = 0;
    std::size_t starP = std::u32string::npos, starT = 0;
    while (t < text.size()) {
        if (p < pattern.size() && pattern[p].kind == GlobItem::Kind::Star) {
            starP = p++;
            starT = t;
        } else if (p < pattern.size() && pattern[p].matches(text[t])) {
            ++p;
            ++t;
        } else if (starP != std::u32string::npos) {
            p = starP + 1;
            t = ++starT;
        } else {
            return false;
        }
    }
    while (p < pattern.size() && pattern[p].kind == GlobItem::Kind::Star) {
        ++p;
    }
    return p == pattern.size();
}

bool hasGlobCharacters(const std::string& text) {
    return text.find_first_of("*?[") != std::string::npos;
}

// ---- the system being selected from ----------------------------------------------------------

enum class Level { Atom, Residue, Segment };

struct Context {
    const Topology& top;
    const Frame* frame;
    std::size_t n;

    std::int32_t residue(std::size_t atom) const { return top.atomResidue[atom]; }
    std::int32_t segment(std::size_t atom) const {
        return top.residueSegment[static_cast<std::size_t>(residue(atom))];
    }

    // the atoms of the residues that pass test(residue index)
    template <class F>
    Mask byResidue(F&& test) const {
        std::vector<std::uint8_t> perResidue(top.residueCount());
        for (std::size_t r = 0; r < perResidue.size(); ++r) {
            perResidue[r] = test(r);
        }
        Mask out(n);
        for (std::size_t a = 0; a < n; ++a) {
            out[a] = perResidue[static_cast<std::size_t>(residue(a))];
        }
        return out;
    }

    // "positions": the current frame, or MDAnalysis' error for a system without coordinates
    const std::vector<float>& positions(const char* missing) const {
        if (!frame) {
            throw AtomSelectionError(missing);
        }
        if (frame->positions.size() != 3 * n) {
            throw std::invalid_argument("frame positions do not match the topology");
        }
        return frame->positions;
    }
};

Mask atomsOfResidues(const Context& ctx, const std::vector<std::uint8_t>& residues) {
    Mask out(ctx.n);
    for (std::size_t a = 0; a < ctx.n; ++a) {
        out[a] = residues[static_cast<std::size_t>(ctx.residue(a))];
    }
    return out;
}

// ---- geometry: periodic minimum-image distances (orthorhombic and triclinic cells) ------------

struct Geometry {
    enum class Mode { None, Ortho, Triclinic } mode = Mode::None;
    std::array<float, 3> box{};
    std::array<float, 3> inverse{};
    std::array<float, 9> tri{};  // row vectors a, b, c (lower triangular); zero when invalid
    std::array<double, 3> spacing{};  // distances between opposite faces of the triclinic cell
    std::optional<std::array<float, 6>> dims;  // the box, when periodic

    // Cell vectors from lengths and angles (degrees): a along x, b in the xy plane, c completing a
    // right-handed cell. All zero for non-positive lengths, angles outside (0, 180) or a flat cell.
    static std::array<float, 9> cellVectors(const std::array<float, 6>& dims) {
        std::array<float, 9> rows{};
        const double la = dims[0], lb = dims[1], lc = dims[2];
        if (!(la > 0 && lb > 0 && lc > 0) || !std::isfinite(la * lb * lc)) {
            return rows;
        }
        std::array<double, 3> cosine{}, sine{};
        for (std::size_t i = 0; i < 3; ++i) {
            const double deg = dims[3 + i];
            if (!(deg > 0 && deg < 180)) {
                return rows;
            }
            if (deg == 90) {
                cosine[i] = 0;
                sine[i] = 1;
            } else {
                const double rad = deg * (std::acos(-1.0) / 180);
                cosine[i] = std::cos(rad);
                sine[i] = std::sin(rad);
            }
        }
        // cosine / sine indexed by alpha (b^c), beta (a^c), gamma (a^b)
        const double cx = cosine[1];
        const double cy = (cosine[0] - cosine[1] * cosine[2]) / sine[2];
        const double czSquared = 1 - cx * cx - cy * cy;
        if (!(czSquared > 0)) {
            return rows;
        }
        const double full[9] = {la, 0, 0, lb * cosine[2], lb * sine[2], 0,
                                lc * cx, lc * cy, lc * std::sqrt(czSquared)};
        for (std::size_t i = 0; i < 9; ++i) {
            rows[i] = static_cast<float>(full[i]);
        }
        if (!(rows[0] > 0 && rows[4] > 0 && rows[8] > 0)) {
            rows.fill(0.0f);
        }
        return rows;
    }

    // group.dimensions when periodic: none without a box or with zero lengths
    static Geometry of(const Frame& frame) {
        Geometry g;
        if (!frame.box) {
            return g;
        }
        const auto& dims = *frame.box;
        if (dims[0] == 0.0f && dims[1] == 0.0f && dims[2] == 0.0f) {
            return g;
        }
        g.dims = dims;
        if (dims[3] == 90.0f && dims[4] == 90.0f && dims[5] == 90.0f) {
            g.mode = Mode::Ortho;
            for (std::size_t i = 0; i < 3; ++i) {
                g.box[i] = dims[i];
                g.inverse[i] = static_cast<float>(1.0 / static_cast<double>(dims[i]));
            }
        } else {
            g.mode = Mode::Triclinic;
            g.tri = cellVectors(dims);
            const auto row = [&g](std::size_t r) {
                return std::array<double, 3>{g.tri[3 * r], g.tri[3 * r + 1], g.tri[3 * r + 2]};
            };
            const double volume = static_cast<double>(g.tri[0]) * g.tri[4] * g.tri[8];
            for (std::size_t k = 0; k < 3 && volume > 0; ++k) {
                const auto u = row((k + 1) % 3), v = row((k + 2) % 3);
                const double nx = u[1] * v[2] - u[2] * v[1];
                const double ny = u[2] * v[0] - u[0] * v[2];
                const double nz = u[0] * v[1] - u[1] * v[0];
                g.spacing[k] = volume / std::sqrt(nx * nx + ny * ny + nz * nz);
            }
        }
        return g;
    }

    double distance(const float* ref, const float* conf) const {
        std::array<double, 3> d{};
        for (std::size_t i = 0; i < 3; ++i) {
            d[i] = static_cast<double>(conf[i]) - static_cast<double>(ref[i]);
        }
        if (mode == Mode::Ortho) {
            for (std::size_t i = 0; i < 3; ++i) {
                if (box[i] > 0.0f && std::isfinite(inverse[i])) {
                    d[i] -= static_cast<double>(box[i]) * std::nearbyint(d[i] * inverse[i]);
                }
            }
        } else if (mode == Mode::Triclinic && tri[0] > 0.0f) {
            return triclinicMinimum(d);
        }
        return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    }

private:
    // Fractional coordinates of a separation (rows of the lower-triangular cell).
    std::array<double, 3> fractional(const std::array<double, 3>& d) const {
        const double z = d[2] / tri[8];
        const double y = (d[1] - z * tri[7]) / tri[4];
        const double x = (d[0] - y * tri[3] - z * tri[6]) / tri[0];
        return {x, y, z};
    }

    // Shortest d + i a + j b + k c. Any image no longer than |d| has fractional coordinate k
    // within |d| / spacing[k] of zero, which bounds the integer shifts to search.
    double triclinicMinimum(std::array<double, 3> d) const {
        auto s = fractional(d);
        if (!std::isfinite(s[0]) || !std::isfinite(s[1]) || !std::isfinite(s[2])) {
            return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        }
        for (std::size_t r = 0; r < 3; ++r) {
            const double n = std::nearbyint(s[r]);
            for (std::size_t i = 0; i < 3; ++i) {
                d[i] -= n * tri[3 * r + i];
            }
        }
        s = fractional(d);
        double best = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        const double radius = std::sqrt(best) * (1 + 1e-9);
        std::array<std::int64_t, 3> lo{}, hi{};
        for (std::size_t r = 0; r < 3; ++r) {
            const double reach = radius / spacing[r];
            lo[r] = static_cast<std::int64_t>(std::ceil(-reach - s[r]));
            hi[r] = static_cast<std::int64_t>(std::floor(reach - s[r]));
        }
        for (std::int64_t i = lo[0]; i <= hi[0]; ++i) {
            for (std::int64_t j = lo[1]; j <= hi[1]; ++j) {
                for (std::int64_t k = lo[2]; k <= hi[2]; ++k) {
                    double sum = 0;
                    for (std::size_t c = 0; c < 3; ++c) {
                        const double v = d[c] + static_cast<double>(i) * tri[c] +
                                         static_cast<double>(j) * tri[3 + c] +
                                         static_cast<double>(k) * tri[6 + c];
                        sum += v * v;
                    }
                    best = std::min(best, sum);
                }
            }
        }
        return std::sqrt(best);
    }
};

// Cells of a uniform grid holding reference points: every point within `reach` of a query lies in
// the 27 cells around the query's cell. Fractional coordinates of the box when periodic.
class CellGrid {
public:
    static std::optional<CellGrid> build(const Geometry& g, const std::vector<float>& ref,
                                         double reach) {
        CellGrid grid;
        const std::size_t count = ref.size() / 3;
        if (!(reach > 0.0) || !std::isfinite(reach) || count == 0) {
            return std::nullopt;
        }
        const double cell = reach * (1.0 + 1e-4) + 1e-4;
        if (g.mode == Geometry::Mode::Ortho || g.mode == Geometry::Mode::Triclinic) {
            std::array<double, 9> m{};
            if (g.mode == Geometry::Mode::Ortho) {
                if (!(g.box[0] > FLT_EPSILON && g.box[1] > FLT_EPSILON && g.box[2] > FLT_EPSILON)) {
                    return std::nullopt;
                }
                m = {g.box[0], 0, 0, 0, g.box[1], 0, 0, 0, g.box[2]};
            } else {
                for (std::size_t i = 0; i < 9; ++i) {
                    m[i] = g.tri[i];
                }
                if (!(m[0] > 0 && m[4] > 0 && m[8] > 0)) {
                    return buildOpen(ref, cell, grid);
                }
            }
            grid.m_periodic = true;
            grid.m_box = m;
            const double volume = m[0] * m[4] * m[8];
            const double bc = std::sqrt(std::pow(m[4] * m[8], 2) + std::pow(m[3] * m[8], 2) +
                                        std::pow(m[3] * m[7] - m[4] * m[6], 2));
            const double ca = m[0] * std::sqrt(m[8] * m[8] + m[7] * m[7]);
            const double heights[3] = {volume / bc, volume / ca, m[8]};
            for (std::size_t i = 0; i < 3; ++i) {
                const double cells = std::floor(heights[i] / cell);
                if (!(cells >= 3)) {
                    return std::nullopt;
                }
                grid.m_dims[i] = static_cast<std::int64_t>(std::min(cells, 256.0));
            }
        } else {
            return buildOpen(ref, cell, grid);
        }
        grid.fill(ref);
        return grid;
    }

    // Calls visit(reference index) for the candidates near p until it returns true.
    template <class F>
    bool any(const float* p, F&& visit) const {
        std::array<std::int64_t, 3> c{};
        if (!cellOf(p, c, true)) {
            return false;
        }
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    std::array<std::int64_t, 3> q = {c[0] + dx, c[1] + dy, c[2] + dz};
                    bool inside = true;
                    for (std::size_t i = 0; i < 3; ++i) {
                        if (m_periodic) {
                            q[i] = (q[i] % m_dims[i] + m_dims[i]) % m_dims[i];
                        } else {
                            inside = inside && q[i] >= 0 && q[i] < m_dims[i];
                        }
                    }
                    if (!inside) {
                        continue;
                    }
                    const auto k =
                        static_cast<std::size_t>((q[0] * m_dims[1] + q[1]) * m_dims[2] + q[2]);
                    for (std::size_t s = m_start[k]; s < m_start[k + 1]; ++s) {
                        if (visit(m_items[s])) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

private:
    static std::optional<CellGrid> buildOpen(const std::vector<float>& ref, double cell,
                                             CellGrid& grid) {
        const std::size_t count = ref.size() / 3;
        std::array<double, 3> lo{}, hi{};
        for (std::size_t i = 0; i < 3; ++i) {
            lo[i] = std::numeric_limits<double>::infinity();
            hi[i] = -std::numeric_limits<double>::infinity();
        }
        for (std::size_t a = 0; a < count; ++a) {
            for (std::size_t i = 0; i < 3; ++i) {
                const double v = ref[3 * a + i];
                if (!std::isfinite(v)) {
                    return std::nullopt;
                }
                lo[i] = std::min(lo[i], v);
                hi[i] = std::max(hi[i], v);
            }
        }
        double total = 1;
        for (std::size_t i = 0; i < 3; ++i) {
            grid.m_dims[i] =
                static_cast<std::int64_t>(std::min(std::floor((hi[i] - lo[i]) / cell) + 1, 1024.0));
            total *= static_cast<double>(grid.m_dims[i]);
        }
        if (total > 4.0 * static_cast<double>(count) + 1e6) {
            return std::nullopt;
        }
        grid.m_low = lo;
        grid.m_cell = cell;
        grid.fill(ref);
        return grid;
    }

    bool cellOf(const float* p, std::array<std::int64_t, 3>& c, bool query) const {
        if (m_periodic) {
            const double z = p[2] / m_box[8];
            const double y = (p[1] - z * m_box[7]) / m_box[4];
            const double x = (p[0] - y * m_box[3] - z * m_box[6]) / m_box[0];
            const double s[3] = {x, y, z};
            for (std::size_t i = 0; i < 3; ++i) {
                if (!std::isfinite(s[i])) {
                    return false;
                }
                const double f = s[i] - std::floor(s[i]);
                c[i] = std::min(m_dims[i] - 1,
                                static_cast<std::int64_t>(f * static_cast<double>(m_dims[i])));
            }
            return true;
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const double v = p[i];
            if (!std::isfinite(v)) {
                return false;
            }
            const double k = std::floor((v - m_low[i]) / m_cell);
            if (query && (k < -1 || k > static_cast<double>(m_dims[i]))) {
                return false;
            }
            c[i] = static_cast<std::int64_t>(k);
        }
        return true;
    }

    void fill(const std::vector<float>& ref) {
        const std::size_t count = ref.size() / 3;
        const auto cells = static_cast<std::size_t>(m_dims[0] * m_dims[1] * m_dims[2]);
        std::vector<std::size_t> cellOfPoint(count);
        m_start.assign(cells + 1, 0);
        for (std::size_t a = 0; a < count; ++a) {
            std::array<std::int64_t, 3> c{};
            cellOf(&ref[3 * a], c, false);
            const auto k = static_cast<std::size_t>((c[0] * m_dims[1] + c[1]) * m_dims[2] + c[2]);
            cellOfPoint[a] = k;
            ++m_start[k + 1];
        }
        for (std::size_t k = 0; k < cells; ++k) {
            m_start[k + 1] += m_start[k];
        }
        m_items.resize(count);
        std::vector<std::size_t> fillAt(m_start.begin(), m_start.end() - 1);
        for (std::size_t a = 0; a < count; ++a) {
            m_items[fillAt[cellOfPoint[a]]++] = a;
        }
    }

    bool m_periodic = false;
    std::array<double, 9> m_box{};
    std::array<double, 3> m_low{};
    double m_cell = 1;
    std::array<std::int64_t, 3> m_dims{1, 1, 1};
    std::vector<std::size_t> m_start;
    std::vector<std::size_t> m_items;
};

// Largest minus smallest coordinate along each axis over the points of both sets (float32).
std::array<float, 3> pointExtent(const std::vector<float>& ref, const std::vector<float>& conf) {
    std::array<float, 3> lo{}, hi{};
    lo.fill(std::numeric_limits<float>::infinity());
    hi.fill(-std::numeric_limits<float>::infinity());
    for (const std::vector<float>* set : {&ref, &conf}) {
        for (std::size_t k = 0; k + 2 < set->size(); k += 3) {
            for (std::size_t i = 0; i < 3; ++i) {
                lo[i] = std::min(lo[i], (*set)[k + i]);
                hi[i] = std::max(hi[i], (*set)[k + i]);
            }
        }
    }
    return {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
}

// MDAnalysis picks a grid search over brute force for a pair search between n reference and m
// configuration points when both sets hold at least 10 points and either n * m >= 1e8 or the
// cutoff is at most 0.3 times the system's size on every axis (float32 arithmetic). Only the grid
// search validates its input, so its errors appear only on that route. Results stay brute force.
void checkGridSearch(const Geometry& g, const std::vector<float>& ref,
                     const std::vector<float>& conf, double cutoff) {
    const std::size_t n = ref.size() / 3, m = conf.size() / 3;
    if (n < 10 || m < 10) {
        return;
    }
    const float cutoff32 = toFloat32(cutoff);
    if (static_cast<double>(n) * static_cast<double>(m) < 1e8) {
        std::array<float, 3> size{};
        if (!g.dims) {
            size = pointExtent(ref, conf);
        } else if (g.mode == Geometry::Mode::Ortho) {
            size = {(*g.dims)[0], (*g.dims)[1], (*g.dims)[2]};
        } else {
            const std::array<float, 9> v = Geometry::cellVectors(*g.dims);
            for (std::size_t i = 0; i < 3; ++i) {
                const float a = v[i], b = v[3 + i], c = v[6 + i];
                size[i] = std::max({a, b, c}) - std::min({a, b, c});
            }
        }
        const float factor = 0.3f;
        for (const float s : size) {
            if (cutoff32 > factor * s) {
                return;  // brute force
            }
        }
    }
    constexpr const char* kZeroBox = "Any of the box dimensions cannot be 0";
    constexpr const char* kNegative = "Cutoff must be positive";
    if (!g.dims) {
        // a non-periodic cubic pseudo-box around the points
        const std::array<float, 3> extent = pointExtent(ref, conf);
        const float largest = std::max({extent[0], extent[1], extent[2]});
        const float edge = std::max(largest, 2.0f * cutoff32) + 2.2f * cutoff32;
        if (edge == 0.0f) {
            throw AtomSelectionError(kZeroBox);
        }
        if (cutoff < 0) {
            throw AtomSelectionError(kNegative);
        }
        return;
    }
    const std::array<float, 6>& dims = *g.dims;
    if (dims[0] == 0.0f || dims[1] == 0.0f || dims[2] == 0.0f) {
        throw AtomSelectionError(kZeroBox);
    }
    if (cutoff < 0) {
        throw AtomSelectionError(kNegative);
    }
    // half the smallest distance between opposite faces as the grid estimates it
    std::array<float, 3> diagonal{}, sine{};
    if (g.mode == Geometry::Mode::Ortho) {
        diagonal = {dims[0], dims[1], dims[2]};
        sine = {1.0f, 1.0f, 1.0f};
    } else {
        const std::array<float, 9> v = Geometry::cellVectors(dims);
        diagonal = {v[0], v[4], v[8]};
        const float radiansPerDegree = static_cast<float>(std::acos(-1.0) / 180.0);
        for (std::size_t i = 0; i < 3; ++i) {
            sine[i] = std::sin(dims[3 + i] * radiansPerDegree);
        }
    }
    const float ax = diagonal[0], by = diagonal[1], cz = diagonal[2];
    const float limit = 0.5f * std::min({by * sine[0], cz * sine[0], ax * sine[1], cz * sine[1],
                                         ax * sine[2], by * sine[2]});
    if (cutoff > static_cast<double>(limit)) {
        throw AtomSelectionError("Cutoff " + pythonFloatRepr(cutoff) + " too large for box (max " +
                                 pythonFloatRepr(static_cast<double>(limit)) + ")");
    }
}

// For every configuration point: whether some reference point lies within each cutoff
// (distance <= cutoff), as capped_distance finds pairs.
std::vector<std::array<bool, 2>> withinCutoffs(const Geometry& g, const std::vector<float>& ref,
                                               const std::vector<float>& conf, double outer,
                                               std::optional<double> inner) {
    const std::size_t nr = ref.size() / 3, nc = conf.size() / 3;
    std::vector<std::array<bool, 2>> out(nc, {false, false});
    auto check = [&](std::size_t c, std::size_t r) {
        const double d = g.distance(&ref[3 * r], &conf[3 * c]);
        out[c][0] = out[c][0] || d <= outer;
        if (!inner) {
            return out[c][0];
        }
        out[c][1] = out[c][1] || d <= *inner;
        return out[c][0] && out[c][1];
    };
    const double reach = inner ? std::max(outer, *inner) : outer;
    std::optional<CellGrid> grid;
    if (static_cast<double>(nr) * static_cast<double>(nc) > 4e6) {
        grid = CellGrid::build(g, ref, reach);
    }
    for (std::size_t c = 0; c < nc; ++c) {
        if (grid) {
            grid->any(&conf[3 * c], [&](std::size_t r) { return check(c, r); });
        } else {
            for (std::size_t r = 0; r < nr && !check(c, r); ++r) {
            }
        }
    }
    return out;
}

// ---- selection tree --------------------------------------------------------------------------

struct Node {
    virtual ~Node() = default;
    virtual Mask apply(const Context& ctx) const = 0;
};
using NodePtr = std::unique_ptr<Node>;

struct AllNode : Node {
    Mask apply(const Context& ctx) const override { return Mask(ctx.n, 1); }
};

struct LogicNode : Node {
    bool isAnd;
    NodePtr left, right;
    LogicNode(bool a, NodePtr l, NodePtr r) : isAnd(a), left(std::move(l)), right(std::move(r)) {}
    Mask apply(const Context& ctx) const override {
        Mask l, r;
        if (isAnd) {  // MDAnalysis applies the right operand first
            r = right->apply(ctx);
            l = left->apply(ctx);
        } else {
            l = left->apply(ctx);
            r = right->apply(ctx);
        }
        for (std::size_t a = 0; a < ctx.n; ++a) {
            l[a] = static_cast<std::uint8_t>(isAnd ? (l[a] & r[a]) : (l[a] | r[a]));
        }
        return l;
    }
};

struct NotNode : Node {
    NodePtr sel;
    explicit NotNode(NodePtr s) : sel(std::move(s)) {}
    Mask apply(const Context& ctx) const override {
        Mask m = sel->apply(ctx);
        for (auto& v : m) {
            v = !v;
        }
        return m;
    }
};

struct GlobalNode : Node {
    NodePtr sel;
    explicit GlobalNode(NodePtr s) : sel(std::move(s)) {}
    Mask apply(const Context& ctx) const override { return sel->apply(ctx); }
};

struct ByResNode : Node {
    NodePtr sel;
    explicit ByResNode(NodePtr s) : sel(std::move(s)) {}
    Mask apply(const Context& ctx) const override {
        const Mask m = sel->apply(ctx);
        std::vector<std::uint8_t> residues(ctx.top.residueCount());
        for (std::size_t a = 0; a < ctx.n; ++a) {
            if (m[a]) {
                residues[static_cast<std::size_t>(ctx.residue(a))] = 1;
            }
        }
        return atomsOfResidues(ctx, residues);
    }
};

struct BondedNode : Node {
    NodePtr sel;
    explicit BondedNode(NodePtr s) : sel(std::move(s)) {}
    Mask apply(const Context& ctx) const override {
        const Mask m = sel->apply(ctx);
        if (!ctx.top.bonds) {
            throw AtomSelectionError("This Universe does not contain bonds information");
        }
        Mask out(ctx.n);
        for (const Bond& b : *ctx.top.bonds) {
            const auto i = static_cast<std::size_t>(b[0]), j = static_cast<std::size_t>(b[1]);
            if (m[j]) {
                out[i] = 1;
            }
            if (m[i]) {
                out[j] = 1;
            }
        }
        return out;
    }
};

// ---- string attributes (fnmatch on each distinct value) --------------------------------------

struct StringField {
    const std::vector<std::uint32_t>* ids = nullptr;  // nullptr: the topology lacks it
    Level level = Level::Atom;
};

StringField stringField(const Topology& top, std::string_view field) {
    if (field == "names") {
        return {&top.names, Level::Atom};
    }
    if (field == "types") {
        return {&top.types, Level::Atom};
    }
    auto opt = [](const auto& values) { return values ? &*values : nullptr; };
    if (field == "elements") {
        return {opt(top.elements), Level::Atom};
    }
    if (field == "altLocs") {
        return {opt(top.altLocs), Level::Atom};
    }
    if (field == "record_types") {
        return {opt(top.recordTypes), Level::Atom};
    }
    if (field == "chainIDs") {
        return {opt(top.chainIDs), Level::Atom};
    }
    if (field == "resnames") {
        return {&top.resnames, Level::Residue};
    }
    if (field == "icodes") {
        return {opt(top.icodes), Level::Residue};
    }
    if (field == "segids") {
        return {&top.segids, Level::Segment};
    }
    return {nullptr, Level::Atom};  // moltypes
}

// Mask of the atoms whose entity's string id passes `test(id)` (called once per distinct id).
template <class F>
Mask maskByStringIds(const Context& ctx, const StringField& field, F&& test) {
    std::vector<std::int8_t> known(ctx.top.strings.size(), -1);
    auto pass = [&](std::uint32_t id) {
        if (known[id] < 0) {
            known[id] = static_cast<std::int8_t>(test(id));
        }
        return known[id] == 1;
    };
    const auto& ids = *field.ids;
    std::vector<std::uint8_t> entity(ids.size());
    for (std::size_t e = 0; e < ids.size(); ++e) {
        entity[e] = pass(ids[e]);
    }
    Mask out(ctx.n);
    for (std::size_t a = 0; a < ctx.n; ++a) {
        switch (field.level) {
        case Level::Atom:
            out[a] = entity[a];
            break;
        case Level::Residue:
            out[a] = entity[static_cast<std::size_t>(ctx.residue(a))];
            break;
        case Level::Segment:
            out[a] = entity[static_cast<std::size_t>(ctx.segment(a))];
            break;
        }
    }
    return out;
}

struct StringNode : Node {
    std::string field;
    std::vector<std::string> patterns;
    Mask apply(const Context& ctx) const override {
        const StringField f = stringField(ctx.top, field);
        if (!f.ids) {
            throw AtomSelectionError("'Topology' object has no attribute '" + field + "'");
        }
        const bool glob = std::any_of(patterns.begin(), patterns.end(), hasGlobCharacters);
        if (!glob) {
            std::set<std::uint32_t> wanted;
            for (const std::string& p : patterns) {
                if (const auto id = ctx.top.strings.find(p)) {
                    wanted.insert(*id);
                }
            }
            return maskByStringIds(ctx, f, [&](std::uint32_t id) { return wanted.count(id) > 0; });
        }
        std::vector<std::vector<GlobItem>> compiled;
        for (const std::string& p : patterns) {
            compiled.push_back(compileGlob(decodeUtf8(p)));
        }
        return maskByStringIds(ctx, f, [&](std::uint32_t id) {
            const std::u32string text = decodeUtf8(ctx.top.strings[id]);
            return std::any_of(compiled.begin(), compiled.end(),
                               [&](const std::vector<GlobItem>& g) { return globMatch(g, text); });
        });
    }
};

// Name sets of the macros (protein, backbone, ...): residue names and, optionally, atom names.
struct MacroNode : Node {
    std::span<const std::string_view> resnames;
    std::span<const std::string_view> names;  // empty: every atom of the residues
    Mask apply(const Context& ctx) const override {
        auto in = [&](std::span<const std::string_view> set, std::uint32_t id) {
            return std::find(set.begin(), set.end(), ctx.top.strings[id]) != set.end();
        };
        Mask out = maskByStringIds(ctx, {&ctx.top.resnames, Level::Residue},
                                   [&](std::uint32_t id) { return in(resnames, id); });
        if (!names.empty()) {
            const Mask byName = maskByStringIds(ctx, {&ctx.top.names, Level::Atom},
                                                [&](std::uint32_t id) { return in(names, id); });
            for (std::size_t a = 0; a < ctx.n; ++a) {
                out[a] = static_cast<std::uint8_t>(out[a] & byName[a]);
            }
        }
        return out;
    }
};

// ---- numeric attributes ----------------------------------------------------------------------

std::string missingAttribute(std::string_view singular) {
    return "This Universe does not contain " + std::string(singular) + " information";
}

// Per-atom integer values of an attribute; std::nullopt when absent.
std::optional<std::vector<std::int64_t>> intValues(const Context& ctx, std::string_view name) {
    const Topology& t = ctx.top;
    std::vector<std::int64_t> out(ctx.n);
    if (name == "index" || name == "bynum") {
        for (std::size_t a = 0; a < ctx.n; ++a) {
            out[a] = static_cast<std::int64_t>(a);
        }
    } else if (name == "resindex") {
        for (std::size_t a = 0; a < ctx.n; ++a) {
            out[a] = ctx.residue(a);
        }
    } else if (name == "segindex") {
        for (std::size_t a = 0; a < ctx.n; ++a) {
            out[a] = ctx.segment(a);
        }
    } else if (name == "resid" || name == "resnum") {
        const auto& values = name == "resid" ? t.resids : t.resnums;
        for (std::size_t a = 0; a < ctx.n; ++a) {
            out[a] = values[static_cast<std::size_t>(ctx.residue(a))];
        }
    } else if (name == "id" && t.ids) {
        out = *t.ids;
    } else if (name == "formalcharge" && t.formalCharges) {
        out = *t.formalCharges;
    } else {
        return std::nullopt;
    }
    return out;
}

std::optional<std::vector<double>> floatValues(const Context& ctx, std::string_view name) {
    const Topology& t = ctx.top;
    const std::optional<std::vector<double>>* values = nullptr;
    if (name == "mass") {
        values = &t.masses;
    } else if (name == "charge") {
        values = &t.charges;
    } else if (name == "occupancy") {
        values = &t.occupancies;
    } else if (name == "tempfactor" || name == "bfactor") {
        values = &t.tempFactors;
    }
    if (!values || !*values) {
        return std::nullopt;
    }
    return **values;
}

std::string singularName(std::string_view keyword) {
    if (keyword == "bfactor") {
        return "tempfactor";
    }
    return std::string(keyword);
}

// index, bynum, resnum, ... <values>: each value an exact integer or an inclusive range
struct IntRangeNode : Node {
    std::string keyword;
    std::vector<PyInt> singles;
    std::vector<std::pair<PyInt, PyInt>> ranges;

    Mask apply(const Context& ctx) const override {
        const auto values = intValues(ctx, keyword);
        if (!values) {
            throw AtomSelectionError(missingAttribute(keyword));
        }
        const std::int64_t shift = keyword == "bynum" ? 1 : 0;  // bynum counts from 1
        Mask out(ctx.n);
        for (std::size_t a = 0; a < ctx.n; ++a) {
            const std::int64_t v = (*values)[a] + shift;
            bool hit = std::any_of(singles.begin(), singles.end(),
                                   [v](const PyInt& x) { return equal(v, x); });
            for (std::size_t k = 0; !hit && k < ranges.size(); ++k) {
                hit = atLeast(v, ranges[k].first) && atMost(v, ranges[k].second);
            }
            out[a] = hit;
        }
        return out;
    }
};

// mass, charge, ... <values>: a single value matches within np.isclose's default tolerances, a
// range inclusively and exactly
struct FloatRangeNode : Node {
    std::string keyword;
    std::vector<double> singles;
    std::vector<std::pair<double, double>> ranges;

    Mask apply(const Context& ctx) const override {
        const auto values = floatValues(ctx, keyword);
        if (!values) {
            throw AtomSelectionError(missingAttribute(singularName(keyword)));
        }
        Mask out(ctx.n);
        for (std::size_t a = 0; a < ctx.n; ++a) {
            const double v = (*values)[a];
            bool hit = std::any_of(singles.begin(), singles.end(),
                                   [v](double x) { return isClose(v, x); });
            for (std::size_t k = 0; !hit && k < ranges.size(); ++k) {
                hit = v >= ranges[k].first && v <= ranges[k].second;
            }
            out[a] = hit;
        }
        return out;
    }
};

struct MissingNode : Node {
    std::string message;
    explicit MissingNode(std::string m) : message(std::move(m)) {}
    Mask apply(const Context&) const override { throw AtomSelectionError(message); }
};

// ---- resid <values>: residue numbers with optional insertion codes --------------------------

struct ResidBound {
    PyInt number;
    std::string icode;  // "" when none was given
};

bool sameNumber(const PyInt& a, const PyInt& b) {
    return a.overflow == b.overflow && (a.overflow != 0 || a.value == b.value);
}

bool below(const PyInt& a, const PyInt& b) {
    if (a.overflow != b.overflow) {
        return a.overflow < b.overflow;
    }
    return a.overflow == 0 && a.value < b.value;
}

// Without an insertion-code attribute: plain integer matching, and no code may be given. With
// one: a bare number means the empty code; a range between different numbers takes the codes from
// the lower bound's code on at its number, every residue strictly between, and the codes up to
// the upper bound's code at its number (applied as such to reversed bounds too).
struct ResidNode : Node {
    std::vector<ResidBound> singles;
    std::vector<std::pair<ResidBound, ResidBound>> ranges;

    Mask apply(const Context& ctx) const override {
        const Topology& t = ctx.top;
        if (!t.icodes) {
            const auto coded = [](const ResidBound& b) { return !b.icode.empty(); };
            const bool anyCode =
                std::any_of(singles.begin(), singles.end(), coded) ||
                std::any_of(ranges.begin(), ranges.end(),
                            [&](const auto& r) { return coded(r.first) || coded(r.second); });
            if (anyCode) {
                throw AtomSelectionError(
                    "Selection specified icodes, while the topology doesn't have any.");
            }
            return ctx.byResidue([&](std::size_t r) {
                const std::int64_t v = t.resids[r];
                return std::any_of(singles.begin(), singles.end(),
                                   [v](const ResidBound& b) { return equal(v, b.number); }) ||
                       std::any_of(ranges.begin(), ranges.end(), [v](const auto& b) {
                           return atLeast(v, b.first.number) && atMost(v, b.second.number);
                       });
            });
        }
        for (const auto& [lo, hi] : ranges) {
            // the residue numbers strictly between the bounds must be countable in 64 bits
            if (below(lo.number, hi.number) &&
                (lo.number.overflow != 0 || hi.number.overflow != 0 ||
                 static_cast<std::uint64_t>(hi.number.value) -
                         static_cast<std::uint64_t>(lo.number.value) >
                     static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))) {
                throw AtomSelectionError("Maximum allowed size exceeded");
            }
        }
        return ctx.byResidue([&](std::size_t r) {
            const std::int64_t v = t.resids[r];
            const std::string& code = t.strings[(*t.icodes)[r]];
            for (const ResidBound& b : singles) {
                if (equal(v, b.number) && code == b.icode) {
                    return true;
                }
            }
            for (const auto& [lo, hi] : ranges) {
                if (sameNumber(lo.number, hi.number)) {
                    if (equal(v, lo.number) && code >= lo.icode && code <= hi.icode) {
                        return true;
                    }
                    continue;
                }
                const bool inside = atLeast(v, lo.number) && !equal(v, lo.number) &&
                                    atMost(v, hi.number) && !equal(v, hi.number);
                if (inside || (equal(v, lo.number) && code >= lo.icode) ||
                    (equal(v, hi.number) && code <= hi.icode)) {
                    return true;
                }
            }
            return false;
        });
    }
};

// ---- atom <segid> <resid> <name> -------------------------------------------------------------

struct AtomNode : Node {
    Token segid;
    PyInt resid;
    Token name;
    Mask apply(const Context& ctx) const override {
        Mask out(ctx.n);
        if (!segid || !name) {
            return out;
        }
        const auto nameId = ctx.top.strings.find(*name);
        const auto segId = ctx.top.strings.find(*segid);
        if (!nameId || !segId) {
            return out;
        }
        for (std::size_t a = 0; a < ctx.n; ++a) {
            const auto r = static_cast<std::size_t>(ctx.residue(a));
            out[a] = ctx.top.names[a] == *nameId && equal(ctx.top.resids[r], resid) &&
                     ctx.top.segids[static_cast<std::size_t>(ctx.segment(a))] == *segId;
        }
        return out;
    }
};

// ---- prop [abs] <attribute> <op> <value> -----------------------------------------------------

std::optional<std::string_view> propAttribute(const Token& prop) {
    if (!prop) {
        return std::string_view{};  // the None entry (attribute None)
    }
    for (const auto& [name, attr] : kPropAttributes) {
        if (name == *prop) {
            return attr;
        }
    }
    return std::nullopt;
}

// bonds, angles, ...: topologies with bonds (parm7) have them in MDAnalysis, MDBridge's do not
bool isConnectivity(std::string_view name) {
    return name == "bonds" || name == "angles" || name == "dihedrals" || name == "impropers" ||
           name == "ureybradleys" || name == "cmaps";
}

bool isPropName(const Token& value) { return propAttribute(value).has_value(); }

// prop [abs] <property> <op> <value>: the per-atom value of an attribute (x, y, z: the current
// positions) compared with a number. Single-precision positions compare with the number rounded to
// float32; == is np.isclose. Text attributes fail as Python's comparisons of str and float do.
struct PropNode : Node {
    Token property;
    bool absolute = false;
    std::string op;  // one of < > <= >= == !=
    double value = 0;

    // values[a] against `threshold` for every atom; `closeEnough` decides ==
    template <class T, class Near>
    Mask compare(const std::vector<T>& values, T threshold, Near closeEnough) const {
        Mask out(values.size());
        for (std::size_t a = 0; a < values.size(); ++a) {
            const T v = absolute ? static_cast<T>(std::fabs(values[a])) : values[a];
            bool hit = v != threshold;
            if (op == "<") {
                hit = v < threshold;
            } else if (op == ">") {
                hit = v > threshold;
            } else if (op == "<=") {
                hit = v <= threshold;
            } else if (op == ">=") {
                hit = v >= threshold;
            } else if (op == "==") {
                hit = closeEnough(v);
            }
            out[a] = hit;
        }
        return out;
    }

    Mask apply(const Context& ctx) const override {
        const auto attr = propAttribute(property);
        if (!attr) {
            throw AtomSelectionError("Expected one of " + std::string(kPropNamesRepr));
        }
        const std::string name = show(property);
        const std::string missing =
            "This Universe does not contain " + std::string(*attr) + " information";
        const auto closeTo = [this](double v) { return isClose(v, value); };
        if (name == "x" || name == "y" || name == "z") {
            const std::vector<float>& xyz = ctx.positions(missing.c_str());
            const std::size_t col = name == "x" ? 0 : name == "y" ? 1 : 2;
            std::vector<float> column(ctx.n);
            for (std::size_t a = 0; a < ctx.n; ++a) {
                column[a] = xyz[3 * a + col];
            }
            return compare(column, toFloat32(value),
                           [this](float v) { return isClose32(v, value); });
        }
        if (name == "mass" || name == "charge" || name == "occupancy" || name == "tempfactor") {
            const auto floats = floatValues(ctx, name);
            if (!floats) {
                throw AtomSelectionError(missing);
            }
            return compare(*floats, value, closeTo);
        }
        if (const auto ints = intValues(ctx, name)) {
            std::vector<double> widened(ints->size());
            for (std::size_t a = 0; a < widened.size(); ++a) {
                widened[a] = static_cast<double>((*ints)[a]);
            }
            return compare(widened, value, closeTo);
        }
        if (name == "type_index" || (isConnectivity(name) && ctx.top.bonds)) {
            throw AtomSelectionError("prop " + name + " is not supported by MDBridge");
        }
        const StringField field = stringField(ctx.top, *attr);
        if (!field.ids || name == "moltype") {
            throw AtomSelectionError(missing);
        }
        if (absolute) {
            throw AtomSelectionError("bad operand type for abs(): 'str'");
        }
        if (op == "!=") {
            return Mask(ctx.n, 1);
        }
        if (op == "==") {
            throw AtomSelectionError("unsupported operand type(s) for -: 'str' and 'float'");
        }
        throw AtomSelectionError("'" + op +
                                 "' not supported between instances of 'str' and 'float'");
    }
};

// ---- same <property> as <selection> ----------------------------------------------------------

struct SameNode : Node {
    Token prop;
    NodePtr sel;

    Mask apply(const Context& ctx) const override {
        const Mask res = sel->apply(ctx);
        if (std::find(res.begin(), res.end(), 1) == res.end()) {
            return Mask(ctx.n);
        }
        const std::string p = show(prop);  // None ends the tokens: expect("as") fails first
        const Topology& t = ctx.top;
        Mask out(ctx.n);
        if (p == "fragment") {
            if (!t.bonds) {
                throw AtomSelectionError("AtomGroup.fragments not available; this requires bonds");
            }
            const auto frag = fragmentIndices(t);
            std::set<std::int32_t> wanted;
            for (std::size_t a = 0; a < ctx.n; ++a) {
                if (res[a]) {
                    wanted.insert(frag[a]);
                }
            }
            for (std::size_t a = 0; a < ctx.n; ++a) {
                out[a] = wanted.count(frag[a]) > 0;
            }
            return out;
        }
        if (p == "x" || p == "y" || p == "z") {
            const auto& xyz = ctx.positions("This Universe has no coordinates");
            const std::size_t col = p == "x" ? 0 : p == "y" ? 1 : 2;
            std::vector<float> values;
            for (std::size_t a = 0; a < ctx.n; ++a) {
                if (res[a] && !std::isnan(xyz[3 * a + col])) {
                    values.push_back(xyz[3 * a + col]);
                }
            }
            std::sort(values.begin(), values.end());
            values.erase(std::unique(values.begin(), values.end()), values.end());
            for (std::size_t a = 0; a < ctx.n; ++a) {
                const float x = xyz[3 * a + col];
                if (std::isnan(x)) {
                    continue;
                }
                // isclose(x, v) needs |x - v| <= 1e-8 + 1e-5 |v|: v within a few 1e-5 |x| of x
                const double w = 3e-5 * std::fabs(static_cast<double>(x)) + 1e-6;
                auto it = std::lower_bound(values.begin(), values.end(), static_cast<float>(x - w));
                if (std::isinf(x)) {
                    it = std::lower_bound(values.begin(), values.end(), x);
                }
                for (; it != values.end() && (*it <= x + w || std::isinf(x)); ++it) {
                    if (isClose32(x, *it)) {
                        out[a] = 1;
                        break;
                    }
                    if (std::isinf(x)) {
                        break;
                    }
                }
            }
            return out;
        }
        const std::string singular = singularName(p);
        if (p == "residue" || p == "segment") {
            const auto vals = *intValues(ctx, p == "residue" ? "resindex" : "segindex");
            return sameInts(ctx, res, vals);
        }
        if (const auto ints = intValues(ctx, p); ints && p != "bynum") {
            return sameInts(ctx, res, *ints);
        }
        if (p == "mass" || p == "charge" || p == "occupancy" || p == "tempfactor" ||
            p == "bfactor") {
            const auto floats = floatValues(ctx, p);
            if (!floats) {
                throw AtomSelectionError(missingAttribute(singular));
            }
            std::set<double> wanted;
            for (std::size_t a = 0; a < ctx.n; ++a) {
                if (res[a] && !std::isnan((*floats)[a])) {
                    wanted.insert((*floats)[a]);
                }
            }
            for (std::size_t a = 0; a < ctx.n; ++a) {
                out[a] = wanted.count((*floats)[a]) > 0;
            }
            return out;
        }
        if (p == "type_index" || (isConnectivity(p) && t.bonds)) {
            throw AtomSelectionError("same " + p + " as ... is not supported by MDBridge");
        }
        const auto attr = propAttribute(prop);
        const StringField field = attr ? stringField(t, *attr) : StringField{};
        if (!field.ids || p == "moltype") {
            throw AtomSelectionError(missingAttribute(singular));
        }
        // the string ids of the selected atoms' entities
        std::set<std::uint32_t> wanted;
        for (std::size_t a = 0; a < ctx.n; ++a) {
            if (res[a]) {
                wanted.insert((*field.ids)[entityOf(ctx, field.level, a)]);
            }
        }
        return maskByStringIds(ctx, field, [&](std::uint32_t id) { return wanted.count(id) > 0; });
    }

    static std::size_t entityOf(const Context& ctx, Level level, std::size_t atom) {
        switch (level) {
        case Level::Residue:
            return static_cast<std::size_t>(ctx.residue(atom));
        case Level::Segment:
            return static_cast<std::size_t>(ctx.segment(atom));
        case Level::Atom:
            break;
        }
        return atom;
    }

    static Mask sameInts(const Context& ctx, const Mask& res,
                         const std::vector<std::int64_t>& vals) {
        std::set<std::int64_t> wanted;
        for (std::size_t a = 0; a < ctx.n; ++a) {
            if (res[a]) {
                wanted.insert(vals[a]);
            }
        }
        Mask out(ctx.n);
        for (std::size_t a = 0; a < ctx.n; ++a) {
            out[a] = wanted.count(vals[a]) > 0;
        }
        return out;
    }
};

// ---- distance keywords -----------------------------------------------------------------------

constexpr const char* kNoDimensions = "AtomGroup has no attribute dimensions. ";

std::vector<float> gather(const std::vector<float>& xyz, const Mask& mask, bool selected) {
    std::vector<float> out;
    for (std::size_t a = 0; a < mask.size(); ++a) {
        if ((mask[a] != 0) == selected) {
            out.insert(out.end(), &xyz[3 * a], &xyz[3 * a] + 3);
        }
    }
    return out;
}

// center_of_geometry() of the selected atoms (float64 mean, then float32 as the selections use it)
std::array<float, 3> centerOfGeometry(const std::vector<float>& xyz, const Mask& mask) {
    double sum[3] = {0, 0, 0};
    std::size_t count = 0;
    for (std::size_t a = 0; a < mask.size(); ++a) {
        if (mask[a]) {
            for (std::size_t i = 0; i < 3; ++i) {
                sum[i] += static_cast<double>(xyz[3 * a + i]);
            }
            ++count;
        }
    }
    std::array<float, 3> out{};
    for (std::size_t i = 0; i < 3; ++i) {
        out[i] = toFloat32(sum[i] / static_cast<double>(count));
    }
    return out;
}

// around / isolayer: atoms outside the selection within `outer` (and not within `inner`) of it
struct AroundNode : Node {
    double outer = 0;
    std::optional<double> inner;  // isolayer
    NodePtr sel;
    Mask apply(const Context& ctx) const override {
        const Mask s = sel->apply(ctx);
        const bool anySel = std::find(s.begin(), s.end(), 1) != s.end();
        const bool anySys = std::find(s.begin(), s.end(), 0) != s.end();
        Mask out(ctx.n);
        if (!anySel || !anySys) {
            return out;
        }
        const auto& xyz = ctx.positions(kNoDimensions);
        const Geometry g = Geometry::of(*ctx.frame);
        const std::vector<float> ref = gather(xyz, s, true);
        const std::vector<float> conf = gather(xyz, s, false);
        checkGridSearch(g, ref, conf, outer);
        if (inner) {
            checkGridSearch(g, ref, conf, *inner);
        }
        const auto hits = withinCutoffs(g, ref, conf, outer, inner);
        std::size_t k = 0;
        for (std::size_t a = 0; a < ctx.n; ++a) {
            if (!s[a]) {
                out[a] = hits[k][0] && !(inner && hits[k][1]);
                ++k;
            }
        }
        return out;
    }
};

// sphzone / sphlayer / point: atoms within a cutoff of one point (beyond `inner` for sphlayer)
struct SphereNode : Node {
    double outer = 0;
    std::optional<double> inner;
    NodePtr sel;  // nullptr: point
    std::array<float, 3> point{};
    Mask apply(const Context& ctx) const override {
        Mask out(ctx.n);
        std::array<float, 3> center = point;
        if (sel) {
            const Mask s = sel->apply(ctx);
            if (std::find(s.begin(), s.end(), 1) == s.end()) {
                return out;
            }
            center = centerOfGeometry(ctx.positions(kNoDimensions), s);
        }
        const auto& xyz = ctx.positions(kNoDimensions);
        const Geometry g = Geometry::of(*ctx.frame);
        for (std::size_t a = 0; a < ctx.n; ++a) {
            const double d = g.distance(center.data(), &xyz[3 * a]);
            out[a] = d <= outer && (!inner || d > *inner);
        }
        return out;
    }
};

// cyzone / cylayer / smarts: parsed as MDAnalysis parses them, not evaluated
struct UnsupportedNode : Node {
    std::string keyword;
    NodePtr sel;
    Mask apply(const Context& ctx) const override {
        if (sel) {
            const Mask s = sel->apply(ctx);
            if (std::find(s.begin(), s.end(), 1) == s.end()) {
                return Mask(ctx.n);
            }
            ctx.positions("This Universe has no coordinates");
        }
        throw AtomSelectionError("'" + keyword + "' selections are not supported by MDBridge");
    }
};

// ---- tokens ----------------------------------------------------------------------------------

// The code point at text[at] and its length in bytes; an invalid sequence counts as one byte that
// is not whitespace.
std::pair<char32_t, std::size_t> codePointAt(std::string_view text, std::size_t at) {
    const auto lead = static_cast<unsigned char>(text[at]);
    std::size_t length = 1;
    char32_t code = lead;
    if (lead >= 0xF0 && lead < 0xF8) {
        length = 4;
        code = lead & 0x07u;
    } else if (lead >= 0xE0) {
        length = lead < 0xF0 ? 3 : 1;
        code = lead & 0x0Fu;
    } else if (lead >= 0xC0) {
        length = 2;
        code = lead & 0x1Fu;
    }
    if (length == 1 || at + length > text.size()) {
        return {lead >= 0x80 ? U'�' : code, 1};
    }
    for (std::size_t k = 1; k < length; ++k) {
        const auto next = static_cast<unsigned char>(text[at + k]);
        if ((next & 0xC0u) != 0x80u) {
            return {U'�', 1};
        }
        code = (code << 6) | (next & 0x3Fu);
    }
    return {code, length};
}

// Words separated by Unicode whitespace; every parenthesis is a word of its own.
std::vector<std::string> tokenize(std::string_view text) {
    std::vector<std::string> words;
    std::string word;
    const auto flush = [&] {
        if (!word.empty()) {
            words.push_back(std::move(word));
            word.clear();
        }
    };
    for (std::size_t at = 0; at < text.size();) {
        const auto [code, length] = codePointAt(text, at);
        if (isPythonSpace(code)) {
            flush();
        } else if (code == U'(' || code == U')') {
            flush();
            words.emplace_back(1, static_cast<char>(code));
        } else {
            word.append(text.substr(at, length));
        }
        at += length;
    }
    flush();
    return words;
}

// ---- keyword table ---------------------------------------------------------------------------

enum class Kind {
    Not,
    Global,
    All,
    Macro,
    Aromatic,
    ByRes,
    Bonded,
    Same,
    Around,
    SphZone,
    SphLayer,
    IsoLayer,
    CyZone,
    CyLayer,
    Point,
    Atom,
    Group,
    Prop,
    Smarts,
    Text,  // glob patterns on a string attribute
    Integer,
    Real,
    Resid,
    Aromaticity,
    Chirality,
};

struct Keyword {
    std::string_view word;
    Kind kind;
    std::string_view field{};  // Text: the attribute (plural name)
    std::span<const std::string_view> resnames{};  // Macro
    std::span<const std::string_view> names{};  // Macro: atom names, empty for all atoms
};

const std::vector<Keyword>& keywords() {
    static const std::vector<Keyword> table = {
        {"not", Kind::Not},
        {"global", Kind::Global},
        {"all", Kind::All},
        {"protein", Kind::Macro, {}, kProteinResnames},
        {"nucleic", Kind::Macro, {}, kNucleicResnames},
        {"water", Kind::Macro, {}, kWaterResnames},
        {"backbone", Kind::Macro, {}, kProteinResnames, kBackboneNames},
        {"nucleicbackbone", Kind::Macro, {}, kNucleicResnames, kNucleicBackboneNames},
        {"nucleicbase", Kind::Macro, {}, kNucleicResnames, kNucleicBaseNames},
        {"nucleicsugar", Kind::Macro, {}, kNucleicResnames, kNucleicSugarNames},
        {"aromatic", Kind::Aromatic},
        {"byres", Kind::ByRes},
        {"bonded", Kind::Bonded},
        {"same", Kind::Same},
        {"around", Kind::Around},
        {"sphzone", Kind::SphZone},
        {"sphlayer", Kind::SphLayer},
        {"isolayer", Kind::IsoLayer},
        {"cyzone", Kind::CyZone},
        {"cylayer", Kind::CyLayer},
        {"point", Kind::Point},
        {"atom", Kind::Atom},
        {"group", Kind::Group},
        {"prop", Kind::Prop},
        {"smarts", Kind::Smarts},
        {"name", Kind::Text, "names"},
        {"type", Kind::Text, "types"},
        {"element", Kind::Text, "elements"},
        {"record_type", Kind::Text, "record_types"},
        {"chainID", Kind::Text, "chainIDs"},
        {"chainid", Kind::Text, "chainIDs"},
        {"altLoc", Kind::Text, "altLocs"},
        {"altloc", Kind::Text, "altLocs"},
        {"resname", Kind::Text, "resnames"},
        {"icode", Kind::Text, "icodes"},
        {"moltype", Kind::Text, "moltypes"},
        {"segid", Kind::Text, "segids"},
        {"index", Kind::Integer},
        {"resindex", Kind::Integer},
        {"segindex", Kind::Integer},
        {"id", Kind::Integer},
        {"bynum", Kind::Integer},
        {"formalcharge", Kind::Integer},
        {"nbindex", Kind::Integer},
        {"resnum", Kind::Integer},
        {"molnum", Kind::Integer},
        {"model", Kind::Integer},
        {"radius", Kind::Real},
        {"tempfactor", Kind::Real},
        {"bfactor", Kind::Real},
        {"mass", Kind::Real},
        {"charge", Kind::Real},
        {"occupancy", Kind::Real},
        {"gbscreen", Kind::Real},
        {"solventradius", Kind::Real},
        {"rmin", Kind::Real},
        {"epsilon", Kind::Real},
        {"rmin14", Kind::Real},
        {"epsilon14", Kind::Real},
        {"resid", Kind::Resid},
        {"aromaticity", Kind::Aromaticity},
        {"chirality", Kind::Chirality},
    };
    return table;
}

const Keyword* findKeyword(const Token& token) {
    if (!token) {
        return nullptr;
    }
    for (const Keyword& k : keywords()) {
        if (k.word == *token) {
            return &k;
        }
    }
    return nullptr;
}

// prop's operand swap found an operator without a reversed form: reported as an unknown token
struct KeyError {};

// ---- parser ----------------------------------------------------------------------------------

// Recursive descent over the tokens with one token of look-ahead. The token list ends with one
// END marker (std::nullopt); keywords that read raw arguments may consume it, after which looking
// at or taking another token fails as reading past the end of Python's deque does.
class Parser {
public:
    explicit Parser(const std::vector<std::string>& words) {
        m_tokens.assign(words.begin(), words.end());
        m_tokens.emplace_back(std::nullopt);
    }

    NodePtr parse() {
        NodePtr node = expression();
        const Token rest = front();
        if (rest) {
            throw AtomSelectionError("Unexpected token at end of selection string: '" + *rest +
                                     "'");
        }
        return node;
    }

private:
    // Deeper nesting than MDAnalysis (which stops at Python's recursion limit) still fails, before
    // the native stack runs out.
    static constexpr int kMaxDepth = 1000;

    class Deeper {
    public:
        explicit Deeper(int& depth) : m_depth(depth) {
            if (++m_depth > kMaxDepth) {
                throw AtomSelectionError("maximum recursion depth exceeded");
            }
        }
        ~Deeper() { --m_depth; }
        Deeper(const Deeper&) = delete;
        Deeper& operator=(const Deeper&) = delete;

    private:
        int& m_depth;
    };

    Token front() const {
        if (m_pos >= m_tokens.size()) {
            throw AtomSelectionError("deque index out of range");
        }
        return m_tokens[m_pos];
    }

    Token pop() {
        if (m_pos >= m_tokens.size()) {
            throw AtomSelectionError("pop from an empty deque");
        }
        return m_tokens[m_pos++];
    }

    static bool isOperator(const Token& token) {
        return token && (*token == "and" || *token == "or");
    }

    // what ends a list of values: END, and / or, parentheses and every keyword
    static bool isKeyword(const Token& token) {
        return !token || isOperator(token) || *token == "(" || *token == ")" ||
               findKeyword(token) != nullptr;
    }

    // primary { (and | or) primary }, combined left to right
    NodePtr expression() {
        const int depthBefore = m_depth;
        NodePtr node = primary();
        try {
            while (isOperator(front())) {
                const bool isAnd = *pop() == "and";
                if (++m_depth > kMaxDepth) {
                    throw AtomSelectionError("maximum recursion depth exceeded");
                }
                NodePtr right = primary();
                node = std::make_unique<LogicNode>(isAnd, std::move(node), std::move(right));
            }
        } catch (...) {
            m_depth = depthBefore;
            throw;
        }
        m_depth = depthBefore;
        return node;
    }

    NodePtr primary() {
        const Deeper deeper(m_depth);
        const Token token = pop();
        if (token && *token == "(") {
            NodePtr inner = expression();
            const Token close = front();
            if (!close || *close != ")") {
                throw AtomSelectionError("Unexpected token: '" + show(close) + "' Expected: ')'");
            }
            pop();
            return inner;
        }
        const Keyword* keyword = findKeyword(token);
        if (!keyword) {
            throw AtomSelectionError("Unknown selection token: '" + show(token) + "'");
        }
        // a keyword's own arguments: value problems are reported as "Selection failed"
        try {
            return keywordSelection(*keyword);
        } catch (const AtomSelectionError&) {
            throw;
        } catch (const KeyError&) {
            throw AtomSelectionError("Unknown selection token: '" + std::string(keyword->word) +
                                     "'");
        } catch (const ValueError& e) {
            throw AtomSelectionError("Selection failed: '" + std::string(e.what()) + "'");
        }
    }

    NodePtr keywordSelection(const Keyword& keyword) {
        switch (keyword.kind) {
        case Kind::Not:
            return std::make_unique<NotNode>(primary());
        case Kind::Global:
            return std::make_unique<GlobalNode>(primary());
        case Kind::All:
            return std::make_unique<AllNode>();
        case Kind::Macro: {
            auto node = std::make_unique<MacroNode>();
            node->resnames = keyword.resnames;
            node->names = keyword.names;
            return node;
        }
        case Kind::Aromatic:
            return std::make_unique<MissingNode>(missingAttribute("aromaticity"));
        case Kind::ByRes:
            return std::make_unique<ByResNode>(expression());
        case Kind::Bonded:
            return std::make_unique<BondedNode>(expression());
        case Kind::Same:
            return sameAs();
        case Kind::Around:
        case Kind::IsoLayer: {
            auto node = std::make_unique<AroundNode>();
            if (keyword.kind == Kind::IsoLayer) {
                node->inner = pyFloat(pop());
            }
            node->outer = pyFloat(pop());
            node->sel = expression();
            return node;
        }
        case Kind::SphZone:
        case Kind::SphLayer: {
            auto node = std::make_unique<SphereNode>();
            if (keyword.kind == Kind::SphLayer) {
                node->inner = pyFloat(pop());
            }
            node->outer = pyFloat(pop());
            node->sel = expression();
            return node;
        }
        case Kind::Point: {
            auto node = std::make_unique<SphereNode>();
            for (float& coordinate : node->point) {
                coordinate = toFloat32(pyFloat(pop()));
            }
            node->outer = pyFloat(pop());
            return node;
        }
        case Kind::CyZone:
        case Kind::CyLayer: {
            const int count = keyword.kind == Kind::CyZone ? 3 : 4;
            for (int k = 0; k < count; ++k) {
                pyFloat(pop());
            }
            auto node = std::make_unique<UnsupportedNode>();
            node->keyword = std::string(keyword.word);
            node->sel = expression();
            return node;
        }
        case Kind::Smarts:
            return smarts();
        case Kind::Atom: {
            auto node = std::make_unique<AtomNode>();
            node->segid = pop();
            node->resid = pyInt(pop());
            node->name = pop();
            return node;
        }
        case Kind::Group: {
            const Token name = pop();
            if (name && *name == "updating") {
                throw AtomSelectionError(
                    "The 'updating' keyword is reserved and cannot be used as a selection group "
                    "name.");
            }
            throw ValueError("Failed to find group: " + show(name));  // no groups are passed
        }
        case Kind::Prop:
            return property();
        case Kind::Text: {
            auto node = std::make_unique<StringNode>();
            node->field = std::string(keyword.field);
            node->patterns = requiredValues();
            return node;
        }
        case Kind::Chirality:
            requiredValues();
            return std::make_unique<MissingNode>(missingAttribute("chirality"));
        case Kind::Aromaticity:
            for (const std::string& value : grabNotKeywords()) {
                const std::string lowered = pythonLower(value);
                if (lowered != "true" && lowered != "false") {
                    throw ValueError("'" + value +
                                     "' is an invalid value for boolean selection. Use 'True' or "
                                     "'False'");
                }
            }
            return std::make_unique<MissingNode>(missingAttribute("aromaticity"));
        case Kind::Integer:
            return integers(std::string(keyword.word));
        case Kind::Real:
            return reals(std::string(keyword.word));
        case Kind::Resid:
            return resids();
        }
        throw std::logic_error("unhandled selection keyword");
    }

    // same <property> as <selection>
    NodePtr sameAs() {
        auto node = std::make_unique<SameNode>();
        node->prop = pop();
        const bool known =
            !node->prop || std::find(std::begin(kSameNames), std::end(kSameNames),
                                     *node->prop) != std::end(kSameNames);
        if (!known) {
            throw ValueError("Unknown same property : " + *node->prop +
                             "Choose one of : " + std::string(kSameNamesRepr));
        }
        const Token as = front();
        if (!as || *as != "as") {
            throw AtomSelectionError("Unexpected token: '" + show(as) + "' Expected: 'as'");
        }
        pop();
        node->sel = expression();
        return node;
    }

    // smarts <pattern>: raw tokens up to a keyword, and / or or END; parentheses are part of the
    // pattern as long as they balance
    NodePtr smarts() {
        int open = 0;
        for (;;) {
            const Token t = front();
            if (t && *t == "(") {
                ++open;
            } else if (t && *t == ")") {
                if (open == 0) {
                    break;
                }
                --open;
            } else if (isKeyword(t)) {
                break;
            }
            pop();
        }
        auto node = std::make_unique<UnsupportedNode>();
        node->keyword = "smarts";
        return node;
    }

    // value keywords: at least one value
    std::vector<std::string> requiredValues() {
        std::vector<std::string> values = grabNotKeywords();
        if (values.empty()) {
            throw ValueError("Unexpected token '" + show(front()) + "'");
        }
        return values;
    }

    // The value tokens after a keyword: every token up to the next keyword, each without one
    // leading backslash (the keyword test sees the token as typed, so "\and" is the value "and").
    std::vector<std::string> grabNotKeywords() {
        std::vector<std::string> values;
        while (!isKeyword(front())) {
            std::string value = *pop();
            if (!value.empty() && value.front() == '\\') {
                value.erase(0, 1);
            }
            values.push_back(std::move(value));
        }
        return values;
    }

    // The values of resid and the numeric keywords, joined into range pieces: a separator token
    // (to, :, -) glues the previous piece, itself and the next token; a token starting with "--",
    // "to" or ":", or following a piece that ends in to, : or -, continues the previous piece.
    std::vector<std::string> rangeValues() {
        const std::vector<std::string> tokens = grabNotKeywords();
        if (tokens.empty()) {
            throw ValueError("Unexpected token: '" + show(front()) + "'");
        }
        const auto startsWith = [](const std::string& text, std::string_view prefix) {
            return text.compare(0, prefix.size(), prefix) == 0;
        };
        const auto endsWith = [](const std::string& text, std::string_view suffix) {
            return text.size() >= suffix.size() &&
                   text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
        };
        const auto joined = [](auto first, auto last) {
            std::string out;
            for (auto it = first; it != last; ++it) {
                out += (it == first ? "" : " ") + *it;
            }
            return out;
        };
        std::vector<std::string> pieces;
        for (std::size_t i = 0; i < tokens.size(); ++i) {
            const std::string& token = tokens[i];
            if (token == "to" || token == ":" || token == "-") {
                if (pieces.empty() || i + 1 == tokens.size()) {
                    throw AtomSelectionError(
                        "Invalid expression given: " + joined(pieces.begin(), pieces.end()) +
                        " " + token + " " +
                        joined(tokens.begin() + static_cast<std::ptrdiff_t>(i + 1), tokens.end()));
                }
                pieces.back() += " " + token + " " + tokens[i + 1];
                ++i;
                continue;
            }
            const bool continues =
                !pieces.empty() &&
                (startsWith(token, "--") || startsWith(token, "to") || startsWith(token, ":") ||
                 endsWith(pieces.back(), "to") || endsWith(pieces.back(), ":") ||
                 endsWith(pieces.back(), "-"));
            if (continues) {
                pieces.back() += " " + token;
            } else {
                pieces.push_back(token);
            }
        }
        return pieces;
    }

    // a whole piece as a Python int, std::nullopt when it is not one
    static std::optional<PyInt> wholeInt(const std::string& piece) {
        try {
            return pyInt(piece);
        } catch (const ValueError&) {
            return std::nullopt;
        }
    }

    NodePtr integers(std::string keyword) {
        static const Rx range = sequence({group(0, intPattern()), rangeSeparator(),
                                          group(1, intPattern())});
        auto node = std::make_unique<IntRangeNode>();
        node->keyword = std::move(keyword);
        for (const std::string& piece : rangeValues()) {
            if (const auto single = wholeInt(piece)) {
                node->singles.push_back(*single);
            } else if (const auto bounds = reMatch(range, 2, piece)) {
                node->ranges.emplace_back(pyInt((*bounds)[0]), pyInt((*bounds)[1]));
            } else {
                throw ValueError("Failed to parse number: " + piece);
            }
        }
        return node;
    }

    NodePtr reals(std::string keyword) {
        static const Rx range = sequence({group(0, floatPattern()), rangeSeparator(),
                                          group(1, floatPattern())});
        auto node = std::make_unique<FloatRangeNode>();
        node->keyword = std::move(keyword);
        for (const std::string& piece : rangeValues()) {
            if (const auto single = tryPythonFloat(piece)) {
                node->singles.push_back(*single);
            } else if (const auto bounds = reMatch(range, 2, piece)) {
                const double lower = pythonFloat((*bounds)[0]);
                node->ranges.emplace_back(lower, pythonFloat((*bounds)[1]));
            } else {
                throw ValueError("Failed to parse number: " + piece);
            }
        }
        return node;
    }

    NodePtr resids() {
        static const Rx single = sequence({group(0, intPattern()), group(1, wordOptional()), end()});
        static const Rx range =
            sequence({group(0, intPattern()), group(1, wordOptional()), rangeSeparator(),
                      group(2, intPattern()), group(3, wordOptional())});
        auto node = std::make_unique<ResidNode>();
        for (const std::string& piece : rangeValues()) {
            if (const auto one = reMatch(single, 2, piece)) {
                node->singles.push_back({pyInt((*one)[0]), (*one)[1]});
            } else if (const auto two = reMatch(range, 4, piece)) {
                node->ranges.emplace_back(ResidBound{pyInt((*two)[0]), (*two)[1]},
                                          ResidBound{pyInt((*two)[2]), (*two)[3]});
            } else {
                throw ValueError("Failed to parse value: " + piece);
            }
        }
        return node;
    }

    // Non-overlapping occurrences of `symbol` in `text`.
    static std::size_t occurrences(const std::string& text, std::string_view symbol) {
        std::size_t count = 0;
        for (std::size_t at = text.find(symbol); at != std::string::npos;
             at = text.find(symbol, at + symbol.size())) {
            ++count;
        }
        return count;
    }

    // prop [abs] <property> <op> <value> (or <value> <op> <property>), from raw tokens; the
    // operator may be written into the property or value token ("x<5", "x <5", "x< 5").
    NodePtr property() {
        static constexpr std::array<std::string_view, 6> kSymbols = {"<=", ">=", "==",
                                                                     "!=", "<",  ">"};
        auto node = std::make_unique<PropNode>();
        Token prop = pop();
        if (prop && *prop == "abs") {
            node->absolute = true;
            prop = pop();
        }
        if (!prop) {
            throw AtomSelectionError("'NoneType' object has no attribute 'split'");
        }
        // an operator symbol found exactly once in the property token splits it
        Token opText;
        for (const std::string_view symbol : kSymbols) {
            if (occurrences(*prop, symbol) == 1) {
                const std::size_t at = prop->find(symbol);
                opText = prop->substr(at);
                prop = prop->substr(0, at);
                break;
            }
        }
        if (!opText) {
            opText = pop();
        }
        if (!opText) {
            throw AtomSelectionError("argument of type 'NoneType' is not iterable");
        }
        // the value may follow the first operator symbol inside the operator text
        Token valueText;
        std::string op = *opText;
        for (const std::string_view symbol : kSymbols) {
            const std::size_t count = occurrences(*opText, symbol);
            if (count == 0) {
                continue;
            }
            if (count > 1) {
                throw ValueError("too many values to unpack (expected 2)");
            }
            const std::size_t after = opText->find(symbol) + symbol.size();
            if (after < opText->size()) {
                op = std::string(symbol);
                valueText = opText->substr(after);
            }
            break;
        }
        if (!valueText) {
            valueText = pop();
        }
        // a property name in the value position: swap, negating the operator
        if (isPropName(valueText)) {
            static const std::map<std::string, std::string, std::less<>> kNegated = {
                {"<", ">="}, {">", "<="}, {"<=", ">"}, {">=", "<"}, {"==", "=="}, {"!=", "!="}};
            const auto it = kNegated.find(op);
            if (it == kNegated.end()) {
                throw KeyError{};
            }
            op = it->second;
            std::swap(prop, valueText);
        }
        if (std::find(kSymbols.begin(), kSymbols.end(), op) == kSymbols.end()) {
            throw ValueError("Invalid operator : '" + op +
                             "' Use one of : 'dict_keys(['>', '<', '>=', '<=', '==', '!='])'");
        }
        node->value = pythonFloat(*valueText);
        node->property = std::move(prop);
        node->op = std::move(op);
        return node;
    }

    std::vector<Token> m_tokens;
    std::size_t m_pos = 0;
    int m_depth = 0;
};

}  // namespace

std::vector<std::int32_t> selectAtoms(const Topology& topology, const Frame* current,
                                      std::string_view selection) {
    if (selection.empty()) {
        return {};  // MDAnalysis warns and selects nothing
    }
    try {
        Parser parser(tokenize(selection));
        const NodePtr tree = parser.parse();
        const Context ctx{topology, current, topology.atomCount()};
        const Mask mask = tree->apply(ctx);
        std::vector<std::int32_t> atoms;
        for (std::size_t a = 0; a < mask.size(); ++a) {
            if (mask[a]) {
                atoms.push_back(static_cast<std::int32_t>(a));
            }
        }
        return atoms;
    } catch (const AtomSelectionError&) {
        throw;
    } catch (const ValueError& e) {
        throw AtomSelectionError(e.what());
    }
}

}  // namespace mdbridge::core
