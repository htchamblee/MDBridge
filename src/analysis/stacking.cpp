// pi-pi stacking between nucleobases and aromatic side chains (mdatk/analysis/stacking.py).
//
// Per frame every aromatic ring gets its centroid, the normal of its least-squares plane and its
// planarity (RMS distance of the ring atoms from that plane). Two rings are stacked when both are
// planar (planarity <= planarity_cutoff), their centroids lie within distance_cutoff, the angle
// between their normals is <= angle_cutoff and the smaller lateral offset of one centroid from the
// other ring's normal axis is <= offset_cutoff. Purines use all nine ring atoms.
//
// Compiled without floating-point contraction (see CMakeLists.txt): centroids, products and sums
// are rounded as numpy rounds them.

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/common.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/core/selections.hpp"

namespace mdbridge::analysis {

namespace {

constexpr double kRadiansPerDegree = 3.141592653589793 / 180.0;  // np.radians: x * (pi / 180)

using RingNames = std::vector<const char*>;

const RingNames kPurine = {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"};
const RingNames kPyrimidine = {"N1", "C2", "N3", "C4", "C5", "C6"};
const RingNames kPhenyl = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"};
const RingNames kIndole = {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"};
const RingNames kImidazole = {"CG", "ND1", "CD2", "CE1", "NE2"};

// AROMATIC_RESIDUES
const RingNames* aromaticRing(const std::string& resname) {
    // resname.upper(): ASCII letters only (MDBridge); Python also maps a few non-ASCII letters
    // ("ı" -> "I", "ſ" -> "S") that no residue name of a real topology uses
    std::string upper = resname;
    for (char& c : upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    static const std::map<std::string, const RingNames*, std::less<>> table = {
        {"PHE", &kPhenyl},    {"TYR", &kPhenyl},    {"TRP", &kIndole},    {"HIS", &kImidazole},
        {"HIE", &kImidazole}, {"HID", &kImidazole}, {"HIP", &kImidazole}, {"HSD", &kImidazole},
        {"HSE", &kImidazole}};
    const auto it = table.find(upper);
    return it == table.end() ? nullptr : it->second;
}

// ctx.box(ts) widened; MDATK's spatial helpers apply _usable_box themselves
std::optional<Box> boxOf(const std::optional<std::array<float, 6>>& dims) {
    if (!dims) {
        return std::nullopt;
    }
    Box out{};
    std::copy(dims->begin(), dims->end(), out.begin());
    return out;
}

double dot3(const Vec3& a, const Vec3& b) { return (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]; }

// np.minimum / np.maximum: NaN propagates
double nanMin(double a, double b) { return (std::isnan(a) || std::isnan(b)) ? (std::isnan(a) ? a : b) : std::min(a, b); }
double nanMax(double a, double b) { return (std::isnan(a) || std::isnan(b)) ? (std::isnan(a) ? a : b) : std::max(a, b); }

class StackingAnalysis final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const AnalysisContext& ctx = context();
        m_p = project().analyses.stacking;
        const core::StackingParams& p = m_p;
        const core::AtomGroup sel = ctx.select(p.selection);
        const std::vector<std::string> names = core::groupNames(project());
        const auto known = [&names](const std::string& name) {
            return std::find(names.begin(), names.end(), name) != names.end();
        };
        const std::vector<std::int32_t> inA =
            known(p.group_a) ? ctx.groupAtoms(p.group_a).residues() : std::vector<std::int32_t>{};
        const std::vector<std::int32_t> inB =
            known(p.group_b) ? ctx.groupAtoms(p.group_b).residues() : std::vector<std::int32_t>{};
        const core::Topology& topology = ctx.topology();
        m_rings.clear();
        m_ringResidues.clear();
        for (const std::int32_t residue : sel.residues()) {
            if (std::optional<std::vector<std::int32_t>> ring = stackingRingAtoms(topology, residue)) {
                m_rings.push_back(std::move(*ring));
                m_ringResidues.push_back(residue);
            }
        }
        if (m_rings.empty()) {
            ctx.warn("stacking: no aromatic rings found in the selection");
        }
        const std::size_t nRings = m_rings.size();
        const std::unordered_set<std::int32_t> setA(inA.begin(), inA.end()), setB(inB.begin(), inB.end());
        m_side.assign(nRings, -1);
        for (std::size_t r = 0; r < nRings; ++r) {
            if (setA.contains(m_ringResidues[r])) {
                m_side[r] = 0;
            } else if (setB.contains(m_ringResidues[r])) {
                m_side[r] = 1;
            }
        }
        // rings grouped by size, sizes in order of first appearance (a dict's order)
        m_bySize.clear();
        for (std::size_t r = 0; r < nRings; ++r) {
            const std::size_t size = m_rings[r].size();
            auto it = std::find_if(m_bySize.begin(), m_bySize.end(),
                                   [size](const SizeGroup& g) { return g.size == size; });
            if (it == m_bySize.end()) {
                m_bySize.push_back({size, {}});
                it = std::prev(m_bySize.end());
            }
            it->rings.push_back(r);
        }
        m_rows.clear();
        m_cols.clear();
        for (std::size_t r = 0; r < nRings; ++r) {
            if (p.scope != "inter" || m_side[r] == 0) {
                m_rows.push_back(static_cast<std::int64_t>(r));
            }
            if (p.scope != "inter" || m_side[r] == 1) {
                m_cols.push_back(static_cast<std::int64_t>(r));
            }
        }
        m_rowOf.assign(nRings, -1);
        m_colOf.assign(nRings, -1);
        for (std::size_t x = 0; x < m_rows.size(); ++x) {
            m_rowOf[static_cast<std::size_t>(m_rows[x])] = static_cast<std::int64_t>(x);
        }
        for (std::size_t x = 0; x < m_cols.size(); ++x) {
            m_colOf[static_cast<std::size_t>(m_cols[x])] = static_cast<std::int64_t>(x);
        }
        const std::size_t n = frameCount();
        m_count.assign(n, 0);
        m_inter.assign(n, 0);
        m_intra.assign(n, 0);
        m_map.reset(n);
        // np.cos(np.radians(x)); numpy's SIMD cos may differ from libm in the last bit
        m_cosCut = std::cos(p.angle_cutoff * kRadiansPerDegree);
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        const std::size_t nRings = m_rings.size();
        if (nRings < 2) {
            return;  // an empty frame of the map
        }
        const core::StackingParams& p = m_p;
        const std::optional<Box> box = boxOf(context().box(frame));
        std::vector<Vec3> cent(nRings), norm(nRings);
        std::vector<double> planar(nRings);
        std::vector<Vec3> positions;
        for (const SizeGroup& group : m_bySize) {
            positions.clear();
            for (const std::size_t r : group.rings) {
                for (const std::int32_t atom : m_rings[r]) {
                    const auto base = 3 * static_cast<std::size_t>(atom);
                    positions.push_back({static_cast<double>(frame.positions.at(base)),
                                         static_cast<double>(frame.positions.at(base + 1)),
                                         static_cast<double>(frame.positions.at(base + 2))});
                }
            }
            const RingGeometry g = ringGeometry(positions, group.size);
            for (std::size_t m = 0; m < group.rings.size(); ++m) {
                cent[group.rings[m]] = g.centroids[m];
                norm[group.rings[m]] = g.normals[m];
                planar[group.rings[m]] = g.planarity[m];
            }
        }
        const auto flat = [&cent](const std::vector<std::int64_t>& which) {
            std::vector<double> out;
            out.reserve(3 * which.size());
            for (const std::int64_t r : which) {
                const Vec3& c = cent[static_cast<std::size_t>(r)];
                out.insert(out.end(), c.begin(), c.end());
            }
            return out;
        };
        // candidate pairs (ring i, ring j); MDBridge's pairsWithin gives them sorted by (i, j)
        std::vector<std::size_t> pi, pj;
        if (p.scope == "inter") {
            const PairList pairs = pairsWithin(flat(m_rows), flat(m_cols), p.distance_cutoff, box);
            for (std::size_t q = 0; q < pairs.size(); ++q) {
                pi.push_back(static_cast<std::size_t>(m_rows[static_cast<std::size_t>(pairs.i[q])]));
                pj.push_back(static_cast<std::size_t>(m_cols[static_cast<std::size_t>(pairs.j[q])]));
            }
        } else {
            std::vector<double> all;
            all.reserve(3 * nRings);
            for (const Vec3& c : cent) {
                all.insert(all.end(), c.begin(), c.end());
            }
            const PairList pairs = pairsWithin(all, all, p.distance_cutoff, box);
            for (std::size_t q = 0; q < pairs.size(); ++q) {
                const auto i = static_cast<std::size_t>(pairs.i[q]);
                const auto j = static_cast<std::size_t>(pairs.j[q]);
                if (i < j && m_ringResidues[i] != m_ringResidues[j]) {
                    pi.push_back(i);
                    pj.push_back(j);
                }
            }
        }
        std::vector<std::int64_t> rowsOut, colsOut;
        std::int64_t inter = 0;
        for (std::size_t q = 0; q < pi.size(); ++q) {
            const std::size_t i = pi[q], j = pj[q];
            Vec3 d{};
            for (std::size_t c = 0; c < 3; ++c) {
                d[c] = cent[j][c] - cent[i][c];
            }
            const Vec3 v = minimumImageVector(d, box);
            const double cosang = std::abs(dot3(norm[i], norm[j]));
            const double projI = dot3(v, norm[i]);
            const double projJ = dot3(v, norm[j]);
            const double vv = dot3(v, v);
            const double offI = std::sqrt(nanMax(vv - projI * projI, 0.0));
            const double offJ = std::sqrt(nanMax(vv - projJ * projJ, 0.0));
            const bool ok = cosang >= m_cosCut && nanMin(offI, offJ) <= p.offset_cutoff &&
                            planar[i] <= p.planarity_cutoff && planar[j] <= p.planarity_cutoff;
            if (!ok) {
                continue;
            }
            rowsOut.push_back(m_rowOf[i]);
            colsOut.push_back(m_colOf[j]);
            if ((m_side[i] == 0 && m_side[j] == 1) || (m_side[i] == 1 && m_side[j] == 0)) {
                ++inter;
            }
        }
        const auto count = static_cast<std::int64_t>(rowsOut.size());
        m_count[k] = count;
        m_inter[k] = inter;
        m_intra[k] = count - inter;
        std::vector<double> ones(rowsOut.size(), 1.0);
        m_map.set(k, {std::move(rowsOut), std::move(colsOut), std::move(ones)});
    }

    core::AnalysisResult finalize() override {
        const core::StackingParams& p = m_p;
        const std::vector<std::string> ringLabels = core::residueLabels(context().topology(), m_ringResidues);
        core::Array rowLabels, colLabels;
        for (const std::int64_t r : m_rows) {
            rowLabels.emplace_back(ringLabels[static_cast<std::size_t>(r)]);
        }
        for (const std::int64_t c : m_cols) {
            colLabels.emplace_back(ringLabels[static_cast<std::size_t>(c)]);
        }
        core::AnalysisResult res = newResult({}, core::Table{{"row_labels", std::move(rowLabels)},
                                                             {"col_labels", std::move(colLabels)},
                                                             {"n_rings", static_cast<std::int64_t>(m_rings.size())},
                                                             {"scope", p.scope}});
        for (auto& entry : m_map.csr().arrays("map")) {
            res.arrays.push_back(std::move(entry));
        }
        res.arrays.emplace_back("count", core::NumArray(m_count));
        res.arrays.emplace_back("inter", core::NumArray(m_inter));
        res.arrays.emplace_back("intra", core::NumArray(m_intra));
        res.addSeries("count", "count", "π–π stacked pairs", "count");
        res.addSeries("inter", "inter", p.group_a + "–" + p.group_b + " stacked pairs", "count");
        if (p.scope != "inter") {
            res.addSeries("intra", "intra", "Intra-group stacked pairs", "count");
        }
        res.addPairmap("map", "map", "row_labels", "col_labels", "Stacking map");
        return res;
    }

private:
    struct SizeGroup {
        std::size_t size;
        std::vector<std::size_t> rings;
    };

    core::StackingParams m_p;
    std::vector<std::vector<std::int32_t>> m_rings;  // ring atoms per ring, in ring order
    std::vector<std::int32_t> m_ringResidues;
    std::vector<int> m_side;  // 0: group A, 1: group B, -1: neither
    std::vector<SizeGroup> m_bySize;
    std::vector<std::int64_t> m_rows, m_cols, m_rowOf, m_colOf;
    std::vector<std::int64_t> m_count, m_inter, m_intra;
    FrameCsr m_map;  // stacked ring pairs per frame
    double m_cosCut = 0.0;
};

class StackingDescriptor final : public AnalysisDescriptor {
public:
    StackingDescriptor()
        : AnalysisDescriptor({"stacking",
                              "π–π stacking",
                              1,
                              false,
                              2,
                              {{"count", "π–π stacked pairs"},
                               {"inter", "{a}–{b} stacked pairs"},
                               {"intra", "Intra-group stacked pairs"},
                               {"map", "Stacking map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.stacking.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "stacking");
    }

    ExpectedOutputs expectedOutputs(const core::Project& project) const override {
        ExpectedOutputs out{{"count", "inter"}, {}, {"map"}};
        if (project.analyses.stacking.scope != "inter") {
            out.series.emplace_back("intra");
        }
        return out;
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::StackingParams& s = project.analyses.stacking;
        return {s.selection, s.group_a, s.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<StackingAnalysis>(*this, context, std::move(instance));
    }
};

}  // namespace

std::optional<std::vector<std::int32_t>> stackingRingAtoms(const core::Topology& topology, std::int32_t residue) {
    // {a.name: a.ix for a in residue.atoms}: the last atom of a repeated name wins
    std::unordered_map<std::string, std::int32_t> byName;
    for (const std::int32_t atom : topology.residueAtoms(residue)) {
        byName[topology.name(atom)] = atom;
    }
    const auto take = [&byName](const RingNames& ring) -> std::optional<std::vector<std::int32_t>> {
        std::vector<std::int32_t> out;
        for (const char* name : ring) {
            const auto it = byName.find(name);
            if (it == byName.end()) {
                return std::nullopt;
            }
            out.push_back(it->second);
        }
        return out;
    };
    for (const RingNames* ring : {&kPurine, &kPyrimidine}) {
        if (auto atoms = take(*ring)) {
            return atoms;
        }
    }
    if (const RingNames* ring = aromaticRing(topology.resname(residue))) {
        return take(*ring);
    }
    return std::nullopt;
}

RingGeometry ringGeometry(std::span<const Vec3> positions, std::size_t ringSize) {
    if (ringSize == 0 || positions.size() % ringSize != 0) {
        throw std::invalid_argument("ringGeometry: the positions do not form rings of the given size");
    }
    const std::size_t nRings = positions.size() / ringSize;
    RingGeometry out;
    out.centroids.resize(nRings);
    std::vector<Vec3> centred(positions.begin(), positions.end());
    for (std::size_t m = 0; m < nRings; ++m) {
        // positions.mean(axis=1): accumulated atom by atom, then divided by the ring size
        Vec3 sum{};
        for (std::size_t a = 0; a < ringSize; ++a) {
            for (std::size_t c = 0; c < 3; ++c) {
                sum[c] += positions[m * ringSize + a][c];
            }
        }
        for (std::size_t c = 0; c < 3; ++c) {
            out.centroids[m][c] = sum[c] / static_cast<double>(ringSize);
        }
        for (std::size_t a = 0; a < ringSize; ++a) {
            for (std::size_t c = 0; c < 3; ++c) {
                centred[m * ringSize + a][c] -= out.centroids[m][c];
            }
        }
    }
    const std::vector<core::SvdRows3> svd = core::svdRows3Batch(centred, ringSize);
    const double root = std::sqrt(static_cast<double>(ringSize));
    out.normals.reserve(nRings);
    out.planarity.reserve(nRings);
    for (const core::SvdRows3& s : svd) {
        out.normals.push_back(s.vt[2]);
        out.planarity.push_back(s.s[2] / root);
    }
    return out;
}

const AnalysisDescriptor& stackingAnalysis() {
    static const StackingDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
