// Mg2+ coordination, bridges and occupancy (mdatk/analysis/mg.py).
//
// * Mg-O contacts: ligand atoms (ligand selection without the ions, ligand elements) within the
//   inner-sphere cutoff of an ion; outer-sphere contacts lie between the inner and outer cutoffs.
// * O-Mg-O bridges: pairs of ligand atoms in one ion's shell (the inner sphere, or everything
//   within the outer cutoff); with `geometry` only pairs whose L-Mg-L angle is near 90 or 180
//   degrees. Bridges between group A and group B, and the ions forming them, are also counted.
// * Occupancy: per residue of the occupancy selection and frame, 0 (no ion within the outer
//   cutoff of its occupancy atoms), 1 (outer sphere) or 2 (inner sphere) from the closest atom.
//
// Compiled without floating-point contraction (see CMakeLists.txt): the cosines are rounded as
// numpy rounds them.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/common.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"

namespace mdbridge::analysis {

namespace {

constexpr double kRadiansPerDegree = 3.141592653589793 / 180.0;  // np.radians: x * (pi / 180)

// atoms.positions.astype(np.float64), x, y, z interleaved
std::vector<double> positionsOf(const core::AtomGroup& atoms, const core::Frame& frame) {
    std::vector<double> out;
    out.reserve(3 * atoms.size());
    for (const std::int32_t atom : atoms.indices()) {
        const auto base = 3 * static_cast<std::size_t>(atom);
        for (std::size_t c = 0; c < 3; ++c) {
            out.push_back(static_cast<double>(frame.positions.at(base + c)));
        }
    }
    return out;
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

// np.isin(atoms.ix, members)
std::vector<bool> memberMask(const core::AtomGroup& atoms, const std::vector<std::int32_t>& members) {
    const std::unordered_set<std::int32_t> set(members.begin(), members.end());
    std::vector<bool> out;
    out.reserve(atoms.size());
    for (const std::int32_t atom : atoms.indices()) {
        out.push_back(set.contains(atom));
    }
    return out;
}

const char* const kOccupancyClasses[] = {"none", "outer sphere", "inner sphere"};

class MgAnalysis final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const AnalysisContext& ctx = context();
        m_p = project().analyses.mg;
        const core::MgParams& p = m_p;
        m_mg = ctx.select(p.mg_selection, true);
        if (m_mg.empty()) {
            ctx.warn("Mg: selection " + core::pythonRepr(core::Value(p.mg_selection)) +
                     " matched no ions; Mg series will be zero");
        }
        const core::AtomGroup ligSel = without(ctx.select(p.ligand_selection), m_mg);
        m_lig = core::elementFilter(ligSel, core::splitCsv(p.ligand_elements));
        const std::vector<std::string> names = core::groupNames(project());
        const auto known = [&names](const std::string& name) {
            return std::find(names.begin(), names.end(), name) != names.end();
        };
        const std::vector<std::int32_t> inA =
            known(p.group_a) ? ctx.groupAtoms(p.group_a).indices() : std::vector<std::int32_t>{};
        const std::vector<std::int32_t> inB =
            known(p.group_b) ? ctx.groupAtoms(p.group_b).indices() : std::vector<std::int32_t>{};
        m_ligA = memberMask(m_lig, inA);
        m_ligB = memberMask(m_lig, inB);
        for (std::size_t k = 0; k < m_ligB.size(); ++k) {
            m_ligB[k] = m_ligB[k] && !m_ligA[k];
        }

        // `p.occupancy_selection or p.ligand_selection`: only the empty string is false in Python
        const std::string& occText = p.occupancy_selection.empty() ? p.ligand_selection : p.occupancy_selection;
        const core::AtomGroup occSel = without(ctx.select(occText), m_mg);
        m_occResidues = residueRows(occSel).residues;  // np.unique(occ_sel.resindices)
        m_occAtoms = core::elementFilter(occSel, core::splitCsv(p.occupancy_elements));
        // np.searchsorted(uniq, occ_atoms.resindices)
        const core::Topology& topology = ctx.topology();
        m_occRow.clear();
        m_occRow.reserve(m_occAtoms.size());
        for (const std::int32_t atom : m_occAtoms.indices()) {
            const std::int32_t residue = topology.atomResidue[static_cast<std::size_t>(atom)];
            const auto it = std::lower_bound(m_occResidues.begin(), m_occResidues.end(), residue);
            m_occRow.push_back(static_cast<std::size_t>(it - m_occResidues.begin()));
        }

        const std::size_t n = frameCount();
        m_contacts.assign(n, 0);
        m_outer.assign(n, 0);
        m_bridges.assign(n, 0);
        m_inter.assign(n, 0);
        m_bridgingIons.assign(n, 0);
        m_occupancy.assign(m_occResidues.size() * n, 0);
        // np.cos(np.radians(x)); numpy's SIMD cos may differ from libm in the last bit (the test
        // data records both)
        m_cos90 = std::cos((90.0 - p.angle_tolerance) * kRadiansPerDegree);  // |cos| <= this ~ 90 deg
        m_cos180 = std::cos((180.0 - p.angle_tolerance) * kRadiansPerDegree);  // cos <= this ~ 180 deg
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        if (m_mg.empty()) {
            return;
        }
        const core::MgParams& p = m_p;
        const std::optional<Box> box = boxOf(context().box(frame));
        const std::vector<double> mpos = positionsOf(m_mg, frame);
        if (!m_lig.empty()) {
            const std::vector<double> lpos = positionsOf(m_lig, frame);
            const PairList pairs = pairsWithin(mpos, lpos, p.outer_cutoff, box);
            const bool innerShell = p.bridge_shell == "inner";
            std::vector<std::int64_t> bi, bj;  // the shell pairs, in pair order
            std::int64_t inner = 0;
            for (std::size_t q = 0; q < pairs.size(); ++q) {
                const bool isInner = pairs.d[q] <= p.inner_cutoff;
                inner += isInner ? 1 : 0;
                if (isInner || !innerShell) {
                    bi.push_back(pairs.i[q]);
                    bj.push_back(pairs.j[q]);
                }
            }
            m_contacts[k] = inner;
            m_outer[k] = static_cast<std::int64_t>(pairs.size()) - inner;
            if (bi.size() > 1) {
                countBridges(k, bi, bj, mpos, lpos, box);
            }
        }
        if (!m_occAtoms.empty()) {
            const std::vector<double> opos = positionsOf(m_occAtoms, frame);
            const PairList pairs = pairsWithin(mpos, opos, p.outer_cutoff, box);
            if (!pairs.empty()) {
                const std::size_t nOcc = m_occResidues.size();
                std::vector<double> best(nOcc, std::numeric_limits<double>::infinity());
                for (std::size_t q = 0; q < pairs.size(); ++q) {  // np.minimum.at (NaN propagates)
                    double& b = best[m_occRow[static_cast<std::size_t>(pairs.j[q])]];
                    const double d = pairs.d[q];
                    b = (std::isnan(b) || std::isnan(d)) ? std::numeric_limits<double>::quiet_NaN() : std::min(b, d);
                }
                const std::size_t n = frameCount();
                for (std::size_t r = 0; r < nOcc; ++r) {
                    std::int8_t cls = 0;
                    if (best[r] <= p.outer_cutoff) {
                        cls = 1;
                    }
                    if (best[r] <= p.inner_cutoff) {
                        cls = 2;
                    }
                    m_occupancy[r * n + k] = cls;
                }
            }
        }
    }

    core::AnalysisResult finalize() override {
        const core::MgParams& p = m_p;
        const std::vector<std::string> occLabels = core::residueLabels(context().topology(), m_occResidues);
        core::AnalysisResult res =
            newResult({}, core::Table{{"occ_labels", core::Array(occLabels.begin(), occLabels.end())},
                                      {"n_ions", static_cast<std::int64_t>(m_mg.size())},
                                      {"group_a", p.group_a},
                                      {"group_b", p.group_b}});
        const std::size_t n = frameCount();
        std::vector<std::int64_t> innerResidues(n, 0);  // (occupancy == 2).sum(axis=0)
        for (std::size_t r = 0; r < m_occResidues.size(); ++r) {
            for (std::size_t k = 0; k < n; ++k) {
                innerResidues[k] += m_occupancy[r * n + k] == 2 ? 1 : 0;
            }
        }
        res.arrays.emplace_back("contacts", core::NumArray(m_contacts));
        res.arrays.emplace_back("outer", core::NumArray(m_outer));
        res.arrays.emplace_back("bridges", core::NumArray(m_bridges));
        res.arrays.emplace_back("inter", core::NumArray(m_inter));
        res.arrays.emplace_back("bridging_ions", core::NumArray(m_bridgingIons));
        res.arrays.emplace_back("occupancy", core::NumArray(m_occupancy, {m_occResidues.size(), n}));
        res.arrays.emplace_back("inner_residues", core::NumArray(std::move(innerResidues)));
        res.addSeries("contacts", "contacts", "Mg–O inner-sphere contacts", "count");
        res.addSeries("outer_contacts", "outer", "Mg–O outer-sphere contacts", "count");
        res.addSeries("bridges", "bridges", "O–Mg–O bridges", "count");
        res.addSeries("inter_bridges", "inter", p.group_a + "–Mg–" + p.group_b + " bridges", "count");
        res.addSeries("bridging_ions", "bridging_ions", "Bridging Mg ions", "count");
        res.addSeries("inner_residues", "inner_residues", "Residues with inner-sphere Mg", "count");
        core::Array classes;
        for (const char* c : kOccupancyClasses) {
            classes.emplace_back(c);
        }
        res.addMatrix("occupancy", "occupancy", "occ_labels", "Mg occupancy", "class",
                      core::Table{{"classes", std::move(classes)}});
        return res;
    }

private:
    // The bridges of frame k from the shell pairs (ion bi[q], ligand bj[q]). The pairs are grouped
    // by ion with a stable sort (np.argsort(kind="stable")), so each ion keeps its ligands in pair
    // order, and every pair of one ion's ligands (np.triu_indices order) is a candidate. The counts
    // do not depend on the pair order.
    void countBridges(std::size_t k, const std::vector<std::int64_t>& bi, const std::vector<std::int64_t>& bj,
                      const std::vector<double>& mpos, const std::vector<double>& lpos,
                      const std::optional<Box>& box) {
        std::vector<std::size_t> order(bi.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&bi](std::size_t x, std::size_t y) { return bi[x] < bi[y]; });
        std::int64_t nBridges = 0, nInter = 0, nIons = 0;
        std::vector<std::size_t> ligs;
        for (std::size_t start = 0; start < order.size();) {
            std::size_t end = start + 1;
            while (end < order.size() && bi[order[end]] == bi[order[start]]) {
                ++end;
            }
            const auto ion = static_cast<std::size_t>(bi[order[start]]);
            ligs.clear();
            for (std::size_t q = start; q < end; ++q) {
                ligs.push_back(static_cast<std::size_t>(bj[order[q]]));
            }
            start = end;
            if (ligs.size() < 2) {
                continue;
            }
            std::int64_t nCross = 0;
            for (std::size_t a = 0; a < ligs.size(); ++a) {
                for (std::size_t b = a + 1; b < ligs.size(); ++b) {
                    const std::size_t l1 = ligs[a], l2 = ligs[b];
                    if (m_p.geometry && !octahedral(ion, l1, l2, mpos, lpos, box)) {
                        continue;
                    }
                    ++nBridges;
                    if ((m_ligA[l1] && m_ligB[l2]) || (m_ligB[l1] && m_ligA[l2])) {
                        ++nCross;
                    }
                }
            }
            nInter += nCross;
            nIons += nCross > 0 ? 1 : 0;
        }
        m_bridges[k] = nBridges;
        m_inter[k] = nInter;
        m_bridgingIons[k] = nIons;
    }

    // c = cos(L1-Mg-L2) of the minimum-image vectors: |c| <= cos90 or c <= cos180
    bool octahedral(std::size_t ion, std::size_t l1, std::size_t l2, const std::vector<double>& mpos,
                    const std::vector<double>& lpos, const std::optional<Box>& box) const {
        Vec3 v1{}, v2{};
        for (std::size_t c = 0; c < 3; ++c) {
            v1[c] = lpos[3 * l1 + c] - mpos[3 * ion + c];
            v2[c] = lpos[3 * l2 + c] - mpos[3 * ion + c];
        }
        v1 = minimumImageVector(v1, box);
        v2 = minimumImageVector(v2, box);
        // (v1 * v2).sum(1) / np.maximum(norm(v1) * norm(v2), 1e-12)
        const double dot = (v1[0] * v2[0] + v1[1] * v2[1]) + v1[2] * v2[2];
        const double n1 = std::sqrt((v1[0] * v1[0] + v1[1] * v1[1]) + v1[2] * v1[2]);
        const double n2 = std::sqrt((v2[0] * v2[0] + v2[1] * v2[1]) + v2[2] * v2[2]);
        const double product = n1 * n2;
        const double denominator = std::isnan(product) ? product : std::max(product, 1e-12);
        const double c = dot / denominator;
        return std::abs(c) <= m_cos90 || c <= m_cos180;
    }

    core::MgParams m_p;
    core::AtomGroup m_mg, m_lig, m_occAtoms;
    std::vector<bool> m_ligA, m_ligB;
    std::vector<std::int32_t> m_occResidues;
    std::vector<std::size_t> m_occRow;
    std::vector<std::int64_t> m_contacts, m_outer, m_bridges, m_inter, m_bridgingIons;
    std::vector<std::int8_t> m_occupancy;  // occupancy residue x frame, row-major
    double m_cos90 = 0.0;
    double m_cos180 = 0.0;
};

class MgDescriptor final : public AnalysisDescriptor {
public:
    MgDescriptor()
        : AnalysisDescriptor({"mg",
                              "Mg coordination",
                              1,
                              false,
                              2,
                              {{"contacts", "Mg–O inner-sphere contacts"},
                               {"outer_contacts", "Mg–O outer-sphere contacts"},
                               {"bridges", "O–Mg–O bridges"},
                               {"inter_bridges", "{a}–Mg–{b} bridges"},
                               {"bridging_ions", "Bridging Mg ions ({a}–{b})"},
                               {"inner_residues", "Residues with inner-sphere Mg"},
                               {"occupancy", "Mg occupancy"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.mg.enabled; }

    core::Table paramsFor(const core::Project& project) const override { return detail::sectionParams(project, "mg"); }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"contacts", "outer_contacts", "bridges", "inter_bridges", "bridging_ions", "inner_residues"},
                {"occupancy"},
                {}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::MgParams& m = project.analyses.mg;
        return {m.mg_selection, m.ligand_selection, m.occupancy_selection, m.group_a, m.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<MgAnalysis>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& mgAnalysis() {
    static const MgDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
