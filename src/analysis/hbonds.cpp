// Hydrogen bonds between two groups (mdatk/analysis/hbonds.py).
//
// A hydrogen bond is counted for every (H, acceptor) pair where the donor heavy atom D and the
// acceptor A are within distance_cutoff and, when angles are used, the D-H...A angle is at least
// angle_cutoff. Donor hydrogens come from the topology's bonds; without bonds they are assigned from
// distances on the first analysed frame. Topologies without hydrogens fall back to donor-acceptor
// distances only.
//
// Deliberate difference: in the distance assignment of hydrogens, a hydrogen goes to its nearest
// candidate donor with MDATK's strict "<" (the first pair wins a tie). MDATK's pair order is
// scipy's tree order; pairsWithin returns pairs sorted by (hydrogen, donor), so an exact distance
// tie goes to the donor with the lowest index here (deterministic), which may be another donor than
// MDATK's.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/common.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"

namespace mdbridge::analysis {

namespace {

core::Array stringArray(const std::vector<std::string>& items) { return core::Array(items.begin(), items.end()); }

std::set<std::string> upperSet(std::string_view csv) {
    std::set<std::string> out;
    for (std::string item : core::splitCsv(csv)) {
        for (char& c : item) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        out.insert(std::move(item));
    }
    return out;
}

// ts.dimensions widened to float64 (std::nullopt without a cell); pairsWithin and
// minimumImageVector apply MDATK's own usability checks.
std::optional<Box> widenedBox(const std::optional<std::array<float, 6>>& dims) {
    if (!dims) {
        return std::nullopt;
    }
    Box box{};
    for (std::size_t n = 0; n < 6; ++n) {
        box[n] = (*dims)[n];
    }
    return box;
}

// np.linalg.norm of one row and the row-wise dot product: ((x0 * y0 + x1 * y1) + x2 * y2)
double dot(const Vec3& x, const Vec3& y) { return (x[0] * y[0] + x[1] * y[1]) + x[2] * y[2]; }

class HBondsAnalysis final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::HBondParams& p = project().analyses.hbonds;
        const core::AtomGroup a = context().groupAtoms(p.group_a);
        const core::AtomGroup b = without(context().groupAtoms(p.group_b), a);
        if (b.empty()) {
            throw core::ValueError("H-bonds: group " + core::pythonRepr(p.group_b) + " lies entirely inside " +
                                   core::pythonRepr(p.group_a));
        }
        m_groupA = p.group_a;
        m_groupB = p.group_b;
        m_ab = a | b;  // sorted union; positions fetched once per frame
        const std::size_t n = m_ab.size();
        const core::Topology& t = context().topology();
        const std::set<std::int32_t> inA(a.indices().begin(), a.indices().end());
        m_side.resize(n);
        for (std::size_t x = 0; x < n; ++x) {
            m_side[x] = inA.contains(m_ab[x]) ? 0 : 1;
        }
        const std::vector<std::string> el = core::atomElements(m_ab);
        m_isH.resize(n);
        m_hasH = false;
        for (std::size_t x = 0; x < n; ++x) {
            m_isH[x] = el[x] == "H";
            m_hasH = m_hasH || m_isH[x];
        }
        const std::set<std::string> donorEl = upperSet(p.donor_elements);
        const std::set<std::string> accEl = upperSet(p.acceptor_elements);
        m_resA = a.residues();
        m_resB = b.residues();
        // residue row of each atom within its own side
        m_row.resize(n);
        for (std::size_t x = 0; x < n; ++x) {
            const std::vector<std::int32_t>& res = m_side[x] == 0 ? m_resA : m_resB;
            const std::int32_t r = t.atomResidue[static_cast<std::size_t>(m_ab[x])];
            m_row[x] = std::lower_bound(res.begin(), res.end(), r) - res.begin();
        }
        m_nB = static_cast<std::int64_t>(m_resB.size());

        m_candDonor.resize(n);
        m_candAcc.resize(n);
        for (std::size_t x = 0; x < n; ++x) {
            m_candDonor[x] = donorEl.contains(el[x]);
            m_candAcc[x] = accEl.contains(el[x]);
        }
        m_hOf.clear();
        bool bondsKnown = false;
        std::vector<core::Bond> bonds;
        try {
            bonds = m_ab.bonds();
            bondsKnown = !bonds.empty();
        } catch (const std::exception&) {
            bonds.clear();
        }
        if (bondsKnown) {
            std::vector<std::int64_t> heavyNeighbors(n, 0);
            std::map<std::int32_t, std::size_t> local;
            for (std::size_t x = 0; x < n; ++x) {
                local.emplace(m_ab[x], x);
            }
            for (const core::Bond& bond : bonds) {
                const auto f1 = local.find(bond[0]);
                const auto f2 = local.find(bond[1]);
                if (f1 == local.end() || f2 == local.end()) {
                    continue;
                }
                const std::size_t l1 = f1->second;
                const std::size_t l2 = f2->second;
                for (const auto& [x, y] : {std::pair{l1, l2}, std::pair{l2, l1}}) {
                    if (m_isH[y]) {
                        hydrogensOf(x).push_back(y);
                    } else {
                        ++heavyNeighbors[x];
                    }
                }
            }
            // nitrogens carrying H or with three heavy neighbours are not acceptors
            for (const auto& [x, hs] : m_hOf) {
                if (el[x] == "N") {
                    m_candAcc[x] = false;
                }
            }
            for (std::size_t x = 0; x < n; ++x) {
                if (el[x] == "N" && heavyNeighbors[x] >= 3) {
                    m_candAcc[x] = false;
                }
            }
        }
        m_useAngle = p.use_angle && m_hasH;
        if (p.use_angle && !m_hasH) {
            context().warn("H-bonds: topology has no hydrogens; using donor-acceptor distance only");
        }
        m_donorsPending = m_hasH && !bondsKnown;
        if (!m_donorsPending) {
            buildDonorLists();
        }
        m_dcut = p.distance_cutoff;
        m_cosCut = std::cos(p.angle_cutoff * (std::numbers::pi / 180.0));  // np.cos(np.radians(x))
        m_count.assign(frameCount(), 0);
        m_map.reset(frameCount());
    }

    // Hydrogens unknown from bonds are guessed on the first analysed frame.
    void prepare(const core::Frame& first) override {
        if (m_donorsPending) {
            guessHydrogens(positions(first), widenedBox(context().box(first)));
            m_donorsPending = false;
        }
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        const std::vector<Vec3> pos = positions(frame);
        const std::optional<Box> box = widenedBox(context().box(frame));
        std::vector<std::int64_t> rowsA, rowsB;
        for (int s = 0; s < 2; ++s) {
            const std::vector<std::size_t>& don = m_donorIdx[s];
            const std::vector<std::size_t>& acc = m_accIdx[1 - s];
            const PairList pairs = pairsWithin(coordinates(pos, don), coordinates(pos, acc), m_dcut, box);
            for (std::size_t e = 0; e < pairs.size(); ++e) {
                const std::size_t d = don[static_cast<std::size_t>(pairs.i[e])];
                const std::size_t ac = acc[static_cast<std::size_t>(pairs.j[e])];
                std::size_t accepted = 1;
                if (m_useAngle) {
                    accepted = 0;
                    for (auto h = static_cast<std::size_t>(m_hStart[d]); h < static_cast<std::size_t>(m_hStart[d + 1]);
                         ++h) {
                        const std::size_t hyd = m_hFlat[h];
                        const Vec3 v1 = minimumImageVector(difference(pos[d], pos[hyd]), box);
                        const Vec3 v2 = minimumImageVector(difference(pos[ac], pos[hyd]), box);
                        const double cosang =
                            dot(v1, v2) / std::max(std::sqrt(dot(v1, v1)) * std::sqrt(dot(v2, v2)), 1e-12);
                        if (cosang <= m_cosCut) {  // angle >= cutoff
                            ++accepted;
                        }
                    }
                }
                for (std::size_t c = 0; c < accepted; ++c) {
                    if (s == 0) {
                        rowsA.push_back(m_row[d]);
                        rowsB.push_back(m_row[ac]);
                    } else {
                        rowsA.push_back(m_row[ac]);
                        rowsB.push_back(m_row[d]);
                    }
                }
            }
        }
        m_count[k] = static_cast<std::int64_t>(rowsA.size());
        if (rowsA.empty()) {
            return;  // an empty frame of the map
        }
        std::vector<std::int64_t> keys(rowsA.size());
        for (std::size_t e = 0; e < keys.size(); ++e) {
            keys[e] = rowsA[e] * m_nB + rowsB[e];
        }
        std::sort(keys.begin(), keys.end());
        FrameEntries entries;
        for (std::size_t e = 0; e < keys.size();) {
            std::size_t end = e;
            while (end < keys.size() && keys[end] == keys[e]) {
                ++end;
            }
            entries.i.push_back(keys[e] / m_nB);
            entries.j.push_back(keys[e] % m_nB);
            entries.v.push_back(static_cast<double>(end - e));
            e = end;
        }
        m_map.set(k, std::move(entries));
    }

    core::AnalysisResult finalize() override {
        const core::HBondParams& p = project().analyses.hbonds;
        const core::Topology& t = context().topology();
        core::AnalysisResult res = newResult({}, {{"labels_a", stringArray(core::residueLabels(t, m_resA))},
                                                  {"labels_b", stringArray(core::residueLabels(t, m_resB))},
                                                  {"group_a", m_groupA},
                                                  {"group_b", m_groupB},
                                                  {"used_angle", m_useAngle},
                                                  {"distance_cutoff", p.distance_cutoff},
                                                  {"angle_cutoff", p.angle_cutoff}});
        core::NamedArrays arrays = m_map.csr().arrays("map");
        const auto offsets = arrays[0].second.values<std::int64_t>();
        core::NumArray perA = perRowCounts(offsets, arrays[1].second.values<std::int32_t>(),
                                           arrays[3].second.values<float>(), m_resA.size());
        core::NumArray perB = perRowCounts(offsets, arrays[2].second.values<std::int32_t>(),
                                           arrays[3].second.values<float>(), m_resB.size());
        for (auto& entry : arrays) {
            res.setArray(entry.first, std::move(entry.second));
        }
        res.setArray("count", core::NumArray(m_count));
        res.setArray("per_res_a", std::move(perA));
        res.setArray("per_res_b", std::move(perB));
        const std::string& ga = m_groupA;
        const std::string& gb = m_groupB;
        res.addSeries("count", "count", ga + "–" + gb + " H-bonds", "count");
        res.addMatrix("per_residue_a", "per_res_a", "labels_a", ga + " residue H-bonds", "count");
        res.addMatrix("per_residue_b", "per_res_b", "labels_b", gb + " residue H-bonds", "count");
        res.addPairmap("map", "map", "labels_a", "labels_b", ga + "–" + gb + " H-bond map");
        return res;
    }

private:
    // the positions of m_ab in `frame`
    std::vector<Vec3> positions(const core::Frame& frame) const {
        const std::size_t n = m_ab.size();
        std::vector<Vec3> pos(n);
        for (std::size_t x = 0; x < n; ++x) {
            const auto base = static_cast<std::size_t>(m_ab[x]) * 3;
            pos[x] = {frame.positions[base], frame.positions[base + 1], frame.positions[base + 2]};
        }
        return pos;
    }

    static Vec3 difference(const Vec3& x, const Vec3& y) { return {x[0] - y[0], x[1] - y[1], x[2] - y[2]}; }

    static std::vector<double> coordinates(const std::vector<Vec3>& pos, const std::vector<std::size_t>& atoms) {
        std::vector<double> out;
        out.reserve(atoms.size() * 3);
        for (const std::size_t x : atoms) {
            out.insert(out.end(), pos[x].begin(), pos[x].end());
        }
        return out;
    }

    // self.h_of.setdefault(x, []): the hydrogens of x, in the order of first insertion of x
    std::vector<std::size_t>& hydrogensOf(std::size_t x) {
        for (auto& [atom, hs] : m_hOf) {
            if (atom == x) {
                return hs;
            }
        }
        return m_hOf.emplace_back(x, std::vector<std::size_t>{}).second;
    }

    void buildDonorLists() {
        const std::size_t n = m_ab.size();
        std::vector<bool> donors = m_candDonor;
        if (m_hasH) {
            std::vector<bool> has(n, false);
            for (const auto& [x, hs] : m_hOf) {
                has[x] = true;
            }
            for (std::size_t x = 0; x < n; ++x) {
                donors[x] = donors[x] && has[x];
            }
        }
        for (int s = 0; s < 2; ++s) {
            m_donorIdx[s].clear();
            m_accIdx[s].clear();
            for (std::size_t x = 0; x < n; ++x) {
                if (m_side[x] == s && donors[x]) {
                    m_donorIdx[s].push_back(x);
                }
                if (m_side[x] == s && m_candAcc[x]) {
                    m_accIdx[s].push_back(x);
                }
            }
        }
        // CSR of hydrogens per donor (local indices)
        std::vector<std::int64_t> counts(n, 0);
        for (const auto& [x, hs] : m_hOf) {
            counts[x] = static_cast<std::int64_t>(hs.size());
        }
        m_hStart.assign(n + 1, 0);
        for (std::size_t x = 0; x < n; ++x) {
            m_hStart[x + 1] = m_hStart[x] + counts[x];
        }
        m_hFlat.assign(static_cast<std::size_t>(m_hStart[n]), 0);
        for (const auto& [x, hs] : m_hOf) {
            std::copy(hs.begin(), hs.end(), m_hFlat.begin() + m_hStart[x]);
        }
    }

    // Hydrogens assigned from distances on the first analysed frame (topology without bonds): each
    // hydrogen goes to the candidate donor nearest within 1.25 Å (strict "<": the first pair wins a
    // tie; pairs come sorted by hydrogen, then donor, so the lowest donor index wins).
    void guessHydrogens(const std::vector<Vec3>& pos, const std::optional<Box>& box) {
        std::vector<std::size_t> heavyIdx, hIdx;
        for (std::size_t x = 0; x < m_ab.size(); ++x) {
            if (m_candDonor[x]) {
                heavyIdx.push_back(x);
            }
            if (m_isH[x]) {
                hIdx.push_back(x);
            }
        }
        const PairList pairs = pairsWithin(coordinates(pos, hIdx), coordinates(pos, heavyIdx), 1.25, box);
        // best[hi] = (hj, dist), in order of first appearance of hi
        std::vector<std::pair<std::int64_t, std::pair<std::int64_t, double>>> best;
        std::map<std::int64_t, std::size_t> slot;
        for (std::size_t e = 0; e < pairs.size(); ++e) {
            const std::int64_t hi = pairs.i[e];
            const std::int64_t hj = pairs.j[e];
            // np.linalg.norm(pos[h] - pos[donor]): the plain difference, no minimum image
            const Vec3 v = difference(pos[hIdx[static_cast<std::size_t>(hi)]], pos[heavyIdx[static_cast<std::size_t>(hj)]]);
            const double dist = std::sqrt(dot(v, v));
            const auto found = slot.find(hi);
            if (found == slot.end()) {
                slot.emplace(hi, best.size());
                best.emplace_back(hi, std::pair{hj, dist});
            } else if (dist < best[found->second].second.second) {
                best[found->second].second = {hj, dist};
            }
        }
        for (const auto& [hi, choice] : best) {
            hydrogensOf(heavyIdx[static_cast<std::size_t>(choice.first)]).push_back(hIdx[static_cast<std::size_t>(hi)]);
        }
        buildDonorLists();
    }

    std::string m_groupA, m_groupB;
    core::AtomGroup m_ab;
    std::vector<int> m_side;
    std::vector<bool> m_isH;
    bool m_hasH = false;
    std::vector<std::int32_t> m_resA, m_resB;
    std::vector<std::int64_t> m_row;
    std::int64_t m_nB = 0;
    std::vector<bool> m_candDonor, m_candAcc;
    std::vector<std::pair<std::size_t, std::vector<std::size_t>>> m_hOf;  // donor -> hydrogens (dict order)
    bool m_useAngle = false;
    bool m_donorsPending = false;
    std::vector<std::size_t> m_donorIdx[2], m_accIdx[2];
    std::vector<std::int64_t> m_hStart;
    std::vector<std::size_t> m_hFlat;
    double m_dcut = 0.0;
    double m_cosCut = 0.0;
    std::vector<std::int64_t> m_count;
    FrameCsr m_map;  // residue H-bond map per frame
};

class HBondsDescriptor final : public AnalysisDescriptor {
public:
    HBondsDescriptor()
        : AnalysisDescriptor({"hbonds",
                              "Hydrogen bonds",
                              1,
                              false,
                              2,
                              {{"count", "{a}–{b} H-bonds"},
                               {"per_residue_a", "{a} residue H-bonds"},
                               {"per_residue_b", "{b} residue H-bonds"},
                               {"map", "{a}–{b} H-bond map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.hbonds.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "hbonds");
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"count"}, {"per_residue_a", "per_residue_b"}, {"map"}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::HBondParams& h = project.analyses.hbonds;
        return {h.group_a, h.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<HBondsAnalysis>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& hbondsAnalysis() {
    static const HBondsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
