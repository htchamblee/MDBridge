// Port of MDATK's mdatk/render/example.py and of the nucleotide geometry of mdatk/demo.py
// (base_template, nucleotide, _level_frame; heavy atoms only), the owner's code. numpy's radians
// is x * (pi / 180), its sin / cos the C library's; 3 x 3 products are written out (numpy's matmul
// of these rotations only adds exact zeros, so the order of its sums does not matter).

#include "mdbridge/render/example.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>

#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::render {

namespace {

using core::AnalysisResult;
using core::NumArray;
using core::Table;
using core::Value;

constexpr double kDegToRad = std::numbers::pi / 180.0;
constexpr double kTwoPi = 2.0 * std::numbers::pi;  // 2 * np.pi
constexpr std::string_view kSequence = "GCAUGGCAUCGA";

char complement(char b) {
    switch (b) {
        case 'A': return 'U';
        case 'U': return 'A';
        case 'G': return 'C';
        default: return 'G';  // C
    }
}

// --- the demo's nucleotide geometry (mdatk/demo.py) -----------------------------------------------

constexpr double kRise = 3.0;
constexpr double kXPyr = 4.046;
constexpr double kXPur = 1.634;
constexpr double kRingR = 1.39;
const double kTwist = 32.7 * kDegToRad;  // np.radians(32.7)

struct Atom {
    std::string name;
    Vec3d xyz;
    std::string el;
};

std::array<double, 2> polar(const std::array<double, 2>& c, double radius, double deg) {
    const double a = deg * kDegToRad;
    return {c[0] + radius * std::cos(a), c[1] + radius * std::sin(a)};
}

// base_template (heavy atoms) and the C1' position
std::vector<Atom> baseHeavyAtoms(char base, std::array<double, 2>& c1) {
    std::vector<Atom> atoms;
    auto add = [&](const std::string& n, const std::array<double, 2>& xy, const std::string& el) {
        atoms.push_back({n, {xy[0], xy[1], 0.0}, el});
    };
    std::map<std::string, std::array<double, 2>> pos;
    if (base == 'C' || base == 'U') {
        const std::array<double, 2> c{-kXPyr, 0.0};
        const std::vector<std::pair<std::string, double>> ring{{"N3", 0}, {"C4", 60}, {"C5", 120},
                                                               {"C6", 180}, {"N1", 240}, {"C2", 300}};
        for (const auto& [n, d] : ring) {
            pos[n] = polar(c, kRingR, d);
        }
        for (const char* n : {"N1", "C2", "N3", "C4", "C5", "C6"}) {
            add(n, pos[n], std::string(1, n[0]));
        }
        add("O2", polar(c, kRingR + 1.23, 300), "O");
        if (base == 'U') {
            add("O4", polar(c, kRingR + 1.23, 60), "O");
        } else {
            add("N4", polar(c, kRingR + 1.34, 60), "N");
        }
        c1 = polar(c, kRingR + 1.47, 240);
    } else {
        const std::array<double, 2> c{-kXPur, 0.0};
        const std::vector<std::pair<std::string, double>> ring{{"N1", 0}, {"C6", 60}, {"C5", 120},
                                                               {"C4", 180}, {"N3", 240}, {"C2", 300}};
        for (const auto& [n, d] : ring) {
            pos[n] = polar(c, kRingR, d);
        }
        const std::array<double, 2> c5{c[0] + -1.871, c[1] + 1.080};
        const double r5 = 1.1824;
        pos["N7"] = polar(c5, r5, 78);
        pos["C8"] = polar(c5, r5, 150);
        pos["N9"] = polar(c5, r5, 222);
        for (const char* n : {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"}) {
            add(n, pos[n], std::string(1, n[0]));
        }
        if (base == 'A') {
            add("N6", polar(c, kRingR + 1.34, 60), "N");
        } else {
            add("O6", polar(c, kRingR + 1.23, 60), "O");
            add("N2", polar(c, kRingR + 1.34, 300), "N");
        }
        c1 = polar(c5, r5 + 1.47, 222);
    }
    return atoms;
}

// nucleotide(base) without hydrogens: atoms and bonds (indices into the atoms) in MDATK's order
std::vector<Atom> heavyNucleotide(char base, std::vector<std::pair<std::size_t, std::size_t>>& bonds) {
    std::array<double, 2> c1{};
    std::vector<Atom> atoms = baseHeavyAtoms(base, c1);
    const Vec3d c1p{c1[0], c1[1], 0.0};
    atoms.push_back({"C1'", c1p, "C"});
    static const std::vector<std::pair<std::string, Vec3d>> sugar{
        {"O4'", {-0.78, 1.08, 0.35}},  {"C2'", {-0.58, -1.07, 0.85}}, {"O2'", {-0.13, -2.37, 0.55}},
        {"C3'", {-1.93, -0.52, 1.25}}, {"C4'", {-2.08, 0.93, 0.85}},  {"O3'", {-2.83, -1.37, 2.00}},
        {"C5'", {-3.28, 1.93, 0.40}},  {"O5'", {-4.18, 0.93, 0.95}},  {"P", {-4.48, 2.03, -0.35}},
        {"OP1", {-5.53, 2.78, -0.95}}, {"OP2", {-5.18, 0.98, -1.15}}};
    for (const auto& [name, off] : sugar) {
        atoms.push_back({name, {c1p[0] + off[0], c1p[1] + off[1], c1p[2] + off[2]}, std::string(1, name[0])});
    }
    std::map<std::string, std::size_t> index;
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        index[atoms[k].name] = k;
    }
    auto bond = [&](const char* x, const char* y) {
        const auto a = index.find(x), b = index.find(y);
        if (a != index.end() && b != index.end()) {
            bonds.emplace_back(a->second, b->second);
        }
    };
    static const std::vector<std::pair<const char*, const char*>> ringBonds{
        {"N1", "C2"}, {"C2", "N3"}, {"N3", "C4"}, {"C4", "C5"}, {"C5", "C6"}, {"C6", "N1"}, {"C5", "N7"},
        {"N7", "C8"}, {"C8", "N9"}, {"N9", "C4"}, {"C2", "O2"}, {"C4", "O4"}, {"C4", "N4"}, {"C6", "N6"},
        {"C6", "O6"}, {"C2", "N2"}};
    for (const auto& [x, y] : ringBonds) {
        bond(x, y);
    }
    bond(base == 'A' || base == 'G' ? "N9" : "N1", "C1'");
    for (const auto& [x, y] : std::vector<std::pair<const char*, const char*>>{
             {"C1'", "O4'"}, {"C1'", "C2'"}, {"C2'", "O2'"}, {"C2'", "C3'"}, {"C3'", "C4'"}, {"C4'", "O4'"},
             {"C3'", "O3'"}, {"C4'", "C5'"}, {"C5'", "O5'"}, {"O5'", "P"}, {"P", "OP1"}, {"P", "OP2"}}) {
        bond(x, y);
    }
    return atoms;
}

// --- arrays ---------------------------------------------------------------------------------------

std::vector<double> flatten(const std::vector<Vec3d>& v, int repeats) {
    std::vector<double> out;
    out.reserve(v.size() * 3 * static_cast<std::size_t>(repeats));
    for (int f = 0; f < repeats; ++f) {
        for (const Vec3d& p : v) {
            out.insert(out.end(), p.begin(), p.end());
        }
    }
    return out;
}

std::vector<std::int64_t> flattenPairs(const std::vector<std::array<std::int64_t, 2>>& v) {
    std::vector<std::int64_t> out;
    for (const auto& p : v) {
        out.push_back(p[0]);
        out.push_back(p[1]);
    }
    return out;
}

std::vector<std::int64_t> arange(std::size_t n) {
    std::vector<std::int64_t> out(n);
    for (std::size_t k = 0; k < n; ++k) {
        out[k] = static_cast<std::int64_t>(k);
    }
    return out;
}

core::Array stringArray(const std::vector<std::string>& v) {
    core::Array out;
    for (const std::string& s : v) {
        out.emplace_back(s);
    }
    return out;
}

core::Array numberLabels(const std::vector<std::int64_t>& v) {
    core::Array out;
    for (std::int64_t r : v) {
        out.emplace_back(std::to_string(r));
    }
    return out;
}

void setMeta(AnalysisResult& r, const std::string& key, Value v) {
    for (auto& [k, val] : r.meta) {
        if (k == key) {
            val = std::move(v);
            return;
        }
    }
    r.meta.emplace_back(key, std::move(v));
}

void setArray(AnalysisResult& r, const std::string& name, NumArray a) {
    if (r.array(name) != nullptr) {
        r.setArray(name, std::move(a));
    } else {
        r.arrays.emplace_back(name, std::move(a));
    }
}

AnalysisResult exampleResult(const std::string& name) {
    AnalysisResult r;
    r.name = name;
    r.key = "example";
    r.meta = Table{{"example", true}};
    return r;
}

// np.clip(x, 0, None): NaN stays, -0.0 becomes +0.0
double clipZero(double x) { return std::isnan(x) || x > 0.0 ? x : 0.0; }

std::string unitOf(std::string_view base, std::string_view out) {
    if (base == "rmsd") {
        return "Å";
    }
    if (base == "energy") {
        return out == "net_force" ? "kcal/mol/Å" : "kcal/mol";
    }
    if (out == "q_native") {
        return "";
    }
    return "count";
}

// 2 pi t / t[-1] (NaN for a single frame, as in MDATK)
std::vector<double> phaseOf(const std::vector<double>& t) {
    std::vector<double> out(t.size(), 0.0);
    for (std::size_t k = 0; k < t.size(); ++k) {
        out[k] = kTwoPi * t[k] / t.back();
    }
    return out;
}

std::vector<double> seriesValues(const std::string& key, const std::string& unit, const std::vector<double>& t) {
    NumpyGenerator rng = exampleRng(key);
    const double phase = rng.uniform(0.0, kTwoPi);
    const std::vector<double> x = phaseOf(t);
    const double freq = rng.uniform(0.8, 2.0);
    std::vector<double> wave(t.size());
    for (std::size_t k = 0; k < t.size(); ++k) {
        wave[k] = std::sin(x[k] * freq + phase);
    }
    std::vector<double> out(t.size());
    if (unit == "count") {
        const auto base = static_cast<double>(rng.integers(8, 30));
        const std::vector<double> noise = rng.normal(0.0, 1.0, t.size());
        for (std::size_t k = 0; k < t.size(); ++k) {
            out[k] = clipZero(std::nearbyint(base + 0.35 * base * wave[k] + noise[k]));
        }
    } else if (unit == "Å") {
        const std::vector<double> noise = rng.normal(0.0, 0.05, t.size());
        for (std::size_t k = 0; k < t.size(); ++k) {
            out[k] = 1.2 + 0.6 * (1.0 + wave[k]) / 2.0 + noise[k];
        }
    } else if (unit.starts_with("kcal")) {
        const std::vector<double> noise = rng.normal(0.0, 4.0, t.size());
        for (std::size_t k = 0; k < t.size(); ++k) {
            out[k] = -120.0 + 35.0 * wave[k] + noise[k];
        }
    } else {
        for (std::size_t k = 0; k < t.size(); ++k) {
            out[k] = 0.8 + 0.15 * wave[k];
        }
    }
    return out;
}

std::vector<std::string> instanceGroups(const core::Project& project, const std::string& instance) {
    try {
        const analysis::AnalysisDescriptor& cls = analysis::analysisFor(instance);
        std::vector<std::string> g = cls.groups(analysis::instanceProject(project, instance));
        if (!g.empty()) {
            return g;
        }
    } catch (const std::out_of_range&) {
    } catch (const core::ValueError&) {
    }
    std::vector<std::string> out;
    for (std::size_t k = 0; k < project.groups.size() && k < 2; ++k) {
        out.push_back(project.groups[k].name);
    }
    return out;
}

AnalysisResult& resultFor(std::vector<AnalysisResult>& results, const std::string& instance) {
    for (AnalysisResult& r : results) {
        if (r.name == instance) {
            return r;
        }
    }
    results.push_back(exampleResult(instance));
    return results.back();
}

void addKey(std::vector<AnalysisResult>& results, const core::Project& project, const std::vector<double>& t,
            const std::string& key, std::string_view kind) {
    if (key.find('/') == std::string::npos) {
        return;
    }
    const auto [instance, out] = analysis::splitKey(key);
    const std::string base = analysis::parseInstance(instance).base;
    AnalysisResult& res = resultFor(results, instance);
    const std::string label = analysis::keyLabel(key, project);
    const std::string arr = "ex_" + std::to_string(res.arrays.size());
    std::vector<std::string> groups = instanceGroups(project, instance);
    groups.emplace_back("A");
    groups.emplace_back("B");
    const std::size_t nFrames = t.size();
    if (kind == "series") {
        const std::string unit = unitOf(base, out);
        setArray(res, arr, NumArray(seriesValues(key, unit, t)));
        res.addSeries(out, arr, label, unit);
    } else if (kind == "matrices") {
        std::vector<std::int64_t> rows;
        if (out.ends_with("_b")) {
            rows = groupResids(project, groups[1]);
        } else if (out == "occupancy") {
            rows = groupResids(project, groups[0]);
            const std::vector<std::int64_t> b = groupResids(project, groups[1]);
            rows.insert(rows.end(), b.begin(), b.end());
        } else {
            rows = groupResids(project, groups[0]);
        }
        setMeta(res, arr + "_rows", numberLabels(rows));
        const std::vector<double> x = phaseOf(t);
        std::vector<double> wave(rows.size() * nFrames);
        for (std::size_t k = 0; k < rows.size(); ++k) {
            for (std::size_t f = 0; f < nFrames; ++f) {
                wave[k * nFrames + f] = std::sin(static_cast<double>(k) / 2.5 + x[f]);
            }
        }
        if (out == "occupancy") {
            std::vector<std::int8_t> cls(wave.size());
            for (std::size_t k = 0; k < wave.size(); ++k) {  // np.digitize(wave, [0.3, 0.8])
                const double w = wave[k];
                cls[k] = static_cast<std::int8_t>(std::isnan(w) ? 2 : (w >= 0.8 ? 2 : (w >= 0.3 ? 1 : 0)));
            }
            setArray(res, arr, NumArray(std::move(cls), {rows.size(), nFrames}));
            res.addMatrix(out, arr, arr + "_rows", label, "class",
                          Table{{"classes", core::Array{Value("none"), Value("outer sphere"), Value("inner sphere")}}});
        } else {
            std::vector<double> v(wave.size());
            for (std::size_t k = 0; k < wave.size(); ++k) {
                v[k] = std::nearbyint(clipZero(4.0 + 4.0 * wave[k]));
            }
            setArray(res, arr, NumArray(std::move(v), {rows.size(), nFrames}));
            res.addMatrix(out, arr, arr + "_rows", label, "count");
        }
    } else if (kind == "pairmaps") {
        const std::vector<std::int64_t> rows = groupResids(project, groups[0]);
        const std::vector<std::int64_t> cols = groupResids(project, groups[1]);
        const auto n = static_cast<std::int64_t>(rows.size());
        const auto m = static_cast<std::int64_t>(cols.size());
        const std::int64_t nn = std::max<std::int64_t>(n, 1), mm = std::max<std::int64_t>(m, 1);
        core::CsrBuilder b;
        for (std::size_t f = 0; f < nFrames; ++f) {
            std::vector<std::int64_t> is, js;
            for (std::int64_t i = 0; i < n; ++i) {
                if (std::sin(static_cast<double>(i) / 2.0 + static_cast<double>(f) / 8.0) > -0.3) {
                    const std::int64_t j = ((m - 1 - i * m / nn) % mm + mm) % mm;  // Python's // and %
                    is.push_back(i);
                    js.push_back(j);
                }
            }
            const std::vector<double> v(is.size(), 3.0);
            b.add(is, js, v);
        }
        for (auto& [name, a] : b.arrays(arr)) {
            setArray(res, name, std::move(a));
        }
        setMeta(res, arr + "_r", numberLabels(rows));
        setMeta(res, arr + "_c", numberLabels(cols));
        res.addPairmap(out, arr, arr + "_r", arr + "_c", label);
    }
}

std::vector<std::string> shownGroups(const core::Project& project) {
    if (!project.visual.structure_groups.empty()) {
        return project.visual.structure_groups;
    }
    std::vector<std::string> names = core::groupNames(project);
    if (names.size() > 2) {
        names.resize(2);
    }
    if (names.empty()) {
        names.emplace_back("A");
    }
    return names;
}

std::int64_t groupPosition(const core::Project& project, const std::string& name) {
    std::int64_t pos = -1;  // {g.name: k}: the last group of that name
    for (std::size_t k = 0; k < project.groups.size(); ++k) {
        if (project.groups[k].name == name) {
            pos = static_cast<std::int64_t>(k);
        }
    }
    return pos;
}

// ("GCAUGGCAUCGA" * 10)[:n] (RN-25: repeated as far as needed)
std::string sequenceOf(std::size_t n) {
    std::string s;
    for (std::size_t k = 0; k < n; ++k) {
        s += kSequence[k % kSequence.size()];
    }
    return s;
}

AnalysisResult exampleBasepairs(const std::vector<std::int64_t>& chain, int nFrames) {
    std::vector<std::int64_t> a, b;
    for (std::size_t k = 0; k < chain.size(); ++k) {
        (chain[k] == chain[0] ? a : b).push_back(static_cast<std::int64_t>(k));
    }
    std::reverse(b.begin(), b.end());
    const std::size_t m = std::min(a.size(), b.size());
    core::CsrBuilder builder;
    for (int f = 0; f < nFrames; ++f) {
        const std::size_t openEnd = (f / 6) % 2 != 0 ? 2 : 0;
        const std::size_t count = m > openEnd ? m - openEnd : 0;  // RN-21: MDATK sliced a[:m - 2] < 0
        const std::vector<std::int64_t> i(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(count));
        const std::vector<std::int64_t> j(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(count));
        builder.add(i, j, std::vector<double>(count, 1.0));
    }
    AnalysisResult res = exampleResult("basepairs");
    res.arrays = builder.arrays("map");
    core::Array labels;
    for (std::size_t k = 0; k < chain.size(); ++k) {
        labels.emplace_back(std::to_string(k));
    }
    res.meta.emplace_back("map_r", labels);
    res.meta.emplace_back("map_c", labels);
    res.addPairmap("map", "map", "map_r", "map_c", "Base pairs");
    return res;
}

std::optional<AnalysisResult> frameBasepairs(const core::Project& project, const TraceSource& src,
                                             const AtomsSource& atoms, int nFrames) {
    const std::vector<std::string>& names = *src.resname;
    const std::size_t nRes = names.size();
    std::vector<chem::AtomIndex> index(nRes);
    for (std::size_t k = 0; k < atoms.names.size() && k < atoms.row.size(); ++k) {
        const std::int64_t row = atoms.row[k];
        if (row >= 0 && static_cast<std::size_t>(row) < nRes) {
            std::string name = atoms.names[k];
            std::replace(name.begin(), name.end(), '*', '\'');
            index[static_cast<std::size_t>(row)].emplace(name, static_cast<std::int64_t>(k));  // setdefault
        }
    }
    const core::BasePairParams& bp = project.analyses.basepairs;
    std::vector<chem::BasePair> found;
    try {
        const chem::PairDetector detector(names, index, bp.cutoff, bp.wobble);
        found = detector.detect(flatten(atoms.coords, 1));
    } catch (const core::ValueError&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
    std::vector<std::int64_t> i, j;
    std::vector<double> v;
    for (const chem::BasePair& p : found) {
        i.push_back(p.i);
        j.push_back(p.j);
        v.push_back(static_cast<double>(p.kind));
    }
    core::CsrBuilder builder;
    for (int f = 0; f < nFrames; ++f) {
        builder.add(i, j, v);
    }
    AnalysisResult res = exampleResult("basepairs");
    res.arrays = builder.arrays("map");
    core::Array labels;
    for (std::size_t k = 0; k < nRes; ++k) {
        labels.emplace_back(std::to_string(k));
    }
    res.meta.emplace_back("map_r", labels);
    res.meta.emplace_back("map_c", labels);
    res.addPairmap("map", "map", "map_r", "map_c", "Base pairs");
    return res;
}

}  // namespace

NumpyGenerator exampleRng(std::string_view key) { return NumpyGenerator(crc32(key.data(), key.size())); }

std::vector<std::int64_t> groupResids(const core::Project& project, std::string_view name, int fallback) {
    std::vector<std::int64_t> out;
    for (const core::ResidueGroup& g : project.groups) {
        if (g.name != name || core::pythonStrip(g.residues).empty()) {
            continue;
        }
        try {
            for (const core::ResidueRange& r : core::parseResidueRanges(g.residues)) {
                // RN-25: stop at the 400 residues used (MDATK listed every number of the range first)
                for (std::int64_t v = r.first; out.size() < kExampleGroupLimit; ++v) {
                    out.push_back(v);
                    if (v == r.last) {
                        break;
                    }
                }
            }
        } catch (const core::ValueError&) {
            out.clear();
        }
        break;
    }
    if (out.empty()) {
        for (int k = 1; k <= fallback; ++k) {
            out.push_back(k);
        }
    }
    return out;
}

TraceSource schematicTrace(const core::Project& project) {
    const std::vector<std::string> names = shownGroups(project);
    TraceSource src;
    src.resname.emplace();
    for (std::size_t s = 0; s < names.size() && s < 2; ++s) {
        const std::vector<std::int64_t> ids = groupResids(project, names[s]);
        std::string seq = sequenceOf(ids.size());
        if (s != 0) {  // the complementary strand, read 5' -> 3'
            std::string comp;
            for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
                comp += complement(*it);
            }
            seq = comp;
        }
        const std::int64_t gpos = groupPosition(project, names[s]);
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const std::size_t level = s == 0 ? k : ids.size() - 1 - k;
            const double ang = (static_cast<double>(level) * 32.7 + (s == 0 ? 0.0 : 150.0)) * kDegToRad;
            src.coords.push_back({9.0 * std::cos(ang), 9.0 * std::sin(ang), 2.8 * static_cast<double>(level)});
            src.chain.push_back(static_cast<std::int64_t>(s));
            src.groupIndex.push_back(gpos);
            src.resid.push_back(ids[k]);
            src.resname->emplace_back(1, seq[k]);
            if (k != 0) {
                const auto last = static_cast<std::int64_t>(src.coords.size()) - 1;
                src.links.push_back({last - 1, last});
            }
        }
    }
    return src;
}

StructurePreview schematicStructure(const core::Project& project) {
    const std::vector<std::string> names = shownGroups(project);
    std::vector<std::int64_t> idsA = groupResids(project, names[0]);
    if (idsA.size() > 40) {
        idsA.resize(40);
    }
    std::vector<std::int64_t> idsB;
    if (names.size() > 1) {
        idsB = groupResids(project, names[1]);
        if (idsB.size() > idsA.size()) {
            idsB.resize(idsA.size());
        }
    }
    const std::string seq = sequenceOf(idsA.size());
    struct Strand {
        std::int64_t index;
        std::string group;
        std::vector<std::int64_t> ids;
        std::string letters;
    };
    std::vector<Strand> strands{{0, names[0], idsA, seq}};
    if (!idsB.empty()) {
        std::string comp;
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            comp += complement(*it);
        }
        strands.push_back({1, names[1], idsB, comp.substr(0, idsB.size())});
    }
    StructurePreview out;
    TraceSource& tr = out.trace;
    AtomsSource& at = out.atoms.emplace();
    tr.resname.emplace();
    for (const Strand& strand : strands) {
        std::optional<std::int64_t> prevO3;
        const std::int64_t gpos = groupPosition(project, strand.group);
        for (std::size_t m = 0; m < strand.ids.size() && m < strand.letters.size(); ++m) {
            const char base = strand.letters[m];
            const std::size_t level = strand.index == 0 ? m : idsA.size() - 1 - m;
            const double t = static_cast<double>(level) * kTwist;
            const double c = std::cos(t), s = std::sin(t);
            const double shiftZ = static_cast<double>(level) * kRise;
            std::vector<std::pair<std::size_t, std::size_t>> nb;
            const std::vector<Atom> nt = heavyNucleotide(base, nb);
            const auto first = static_cast<std::int64_t>(at.coords.size());
            const auto rr = static_cast<std::int64_t>(tr.coords.size());
            std::map<std::string, std::int64_t> idx;
            for (std::size_t k = 0; k < nt.size(); ++k) {
                Vec3d v = nt[k].xyz;
                if (strand.index != 0) {
                    v = {-v[0], v[1], -v[2]};  // DYAD
                }
                at.coords.push_back({c * v[0] + -s * v[1] + 0.0, s * v[0] + c * v[1] + 0.0, v[2] + shiftZ});
                at.names.push_back(nt[k].name);
                at.elements.push_back(nt[k].el);
                at.resnames.emplace_back(1, base);
                at.row.push_back(rr);
                idx[nt[k].name] = first + static_cast<std::int64_t>(k);
            }
            for (const auto& [i, j] : nb) {
                at.bonds.push_back({first + static_cast<std::int64_t>(i), first + static_cast<std::int64_t>(j)});
            }
            if (prevO3 && idx.contains("P")) {
                at.bonds.push_back({*prevO3, idx["P"]});
            }
            prevO3 = idx.contains("O3'") ? std::optional<std::int64_t>(idx["O3'"]) : std::nullopt;
            tr.coords.push_back(at.coords[static_cast<std::size_t>(idx["C4'"])]);
            tr.resname->emplace_back(1, base);
            tr.chain.push_back(strand.index);
            tr.groupIndex.push_back(gpos);
            tr.resid.push_back(strand.ids[m]);
            if (m != 0) {
                tr.links.push_back({rr - 1, rr});
            }
        }
    }
    return out;
}

core::Dataset exampleDataset(const core::Project& project, const StructurePreview* preview, int nFrames) {
    if (nFrames < 1) {
        throw core::ValueError("the example dataset needs at least one frame");
    }
    const auto n = static_cast<std::size_t>(nFrames);
    const double dt = project.inputs.dt_ps.has_value() && *project.inputs.dt_ps != 0.0 ? *project.inputs.dt_ps : 10.0;
    std::vector<double> t(n);
    for (std::size_t k = 0; k < n; ++k) {
        t[k] = static_cast<double>(k) * dt;
    }
    std::vector<AnalysisResult> results;
    bool structural = false;
    for (const core::PanelSpec& spec : project.layout.panels) {
        structural = structural || spec.type == "structure" || spec.type == "sequence" || spec.type == "secondary";
        const core::PanelType* cls = core::findPanelType(spec.type);
        if (cls == nullptr) {
            continue;
        }
        for (const std::string& key : spec.series) {
            addKey(results, project, t, key, cls->dataKind);
        }
    }
    if (structural) {
        std::optional<StructurePreview> schematic;
        const TraceSource* src = nullptr;
        const AtomsSource* atoms = nullptr;
        if (preview != nullptr) {
            src = &preview->trace;
            atoms = preview->atoms ? &*preview->atoms : nullptr;
        } else if (!project.visual.representations.empty()) {
            schematic = schematicStructure(project);
            src = &schematic->trace;
            atoms = &*schematic->atoms;
        } else {
            schematic = StructurePreview{schematicTrace(project), std::nullopt};
            src = &schematic->trace;
        }
        const std::size_t nRes = src->coords.size();
        AnalysisResult trace = exampleResult("trace");
        trace.meta.emplace_back("schematic", preview == nullptr);
        trace.arrays.emplace_back("coords", NumArray(flatten(src->coords, nFrames), {n, nRes, 3}));
        trace.arrays.emplace_back("ion_coords", NumArray(flatten(src->ions, nFrames), {n, src->ions.size(), 3}));
        trace.arrays.emplace_back("chain", NumArray(src->chain));
        trace.arrays.emplace_back("group_index", NumArray(src->groupIndex));
        const std::size_t nGroups = src->groupMember ? src->groupMember->size() : project.groups.size();
        std::vector<bool> member;
        for (std::size_t g = 0; g < nGroups; ++g) {
            for (std::size_t r = 0; r < nRes; ++r) {
                member.push_back(src->groupMember ? (*src->groupMember)[g][r]
                                                  : src->groupIndex[r] == static_cast<std::int64_t>(g));
            }
        }
        trace.arrays.emplace_back("group_member", NumArray(std::move(member), {nGroups, nRes}));
        trace.arrays.emplace_back("is_point", NumArray(std::vector<bool>(nRes, true)));
        trace.arrays.emplace_back("links", NumArray(flattenPairs(src->links), {src->links.size(), 2}));
        trace.arrays.emplace_back("resindex", NumArray(arange(nRes)));
        trace.arrays.emplace_back("resid", NumArray(src->resid));
        if (src->resname) {
            trace.arrays.emplace_back("resname", NumArray::strings(*src->resname, 8));
            std::vector<std::string> chainName;
            if (src->chainName) {
                chainName = *src->chainName;
            } else {
                for (std::int64_t c : src->chain) {
                    chainName.emplace_back(1, static_cast<char>('A' + ((c % 26) + 26) % 26));
                }
            }
            trace.arrays.emplace_back("chain_name", NumArray(chainName));
        }
        trace.meta.emplace_back("residue_labels", stringArray(src->labels));
        results.push_back(std::move(trace));
        std::optional<AnalysisResult> pairs;
        if (preview != nullptr && atoms != nullptr && src->resname) {
            pairs = frameBasepairs(project, *src, *atoms, nFrames);
        }
        results.push_back(pairs ? std::move(*pairs) : exampleBasepairs(src->chain, nFrames));
        if (atoms != nullptr) {
            AnalysisResult a = exampleResult("atoms");
            a.meta.emplace_back("atom_names", stringArray(atoms->names));
            a.meta.emplace_back("elements", stringArray(atoms->elements));
            a.meta.emplace_back("resnames", stringArray(atoms->resnames));
            std::vector<float> xyz;
            for (double v : flatten(atoms->coords, nFrames)) {
                xyz.push_back(static_cast<float>(v));
            }
            a.arrays.emplace_back("coords", NumArray(std::move(xyz), {n, atoms->coords.size(), 3}));
            a.arrays.emplace_back("row", NumArray(atoms->row));
            a.arrays.emplace_back("bonds", NumArray(flattenPairs(atoms->bonds), {atoms->bonds.size(), 2}));
            a.arrays.emplace_back("resindex", NumArray(arange(nRes)));
            results.push_back(std::move(a));
        }
        // proximity: every residue drifts in and out of contact so glow modes are visible
        std::vector<float> dist(nRes * n);
        for (std::size_t k = 0; k < nRes; ++k) {
            for (std::size_t f = 0; f < n; ++f) {
                dist[k * n + f] = static_cast<float>(
                    5.0 + 4.0 * std::sin(static_cast<double>(k) / 3.0 + static_cast<double>(f) / 9.0));
            }
        }
        AnalysisResult prox = exampleResult("proximity");
        prox.arrays.emplace_back("min_distance", NumArray(std::move(dist), {nRes, n}));
        prox.arrays.emplace_back("member", NumArray(std::vector<bool>(nRes, false)));
        prox.arrays.emplace_back("resindex", NumArray(arange(nRes)));
        results.push_back(std::move(prox));
    }
    std::vector<std::int64_t> frames = arange(n);
    return core::Dataset(std::move(frames), std::move(t), std::move(results), project.inputs.time_unit);
}

}  // namespace mdbridge::render
