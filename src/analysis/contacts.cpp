// Heavy-atom and residue contacts between two groups (mdatk/analysis/contacts.py).

#include <algorithm>
#include <cstdint>
#include <memory>
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

// positions of `atoms` (x, y, z per atom, float32 as the frame holds them)
std::vector<float> gather(const core::Frame& frame, const std::vector<std::int32_t>& atoms) {
    std::vector<float> out;
    out.reserve(atoms.size() * 3);
    for (const std::int32_t atom : atoms) {
        const auto base = static_cast<std::size_t>(atom) * 3;
        out.insert(out.end(), frame.positions.begin() + static_cast<std::ptrdiff_t>(base),
                   frame.positions.begin() + static_cast<std::ptrdiff_t>(base + 3));
    }
    return out;
}

core::Array stringArray(const std::vector<std::string>& items) { return core::Array(items.begin(), items.end()); }

class ContactsAnalysis final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::ContactParams& p = project().analyses.contacts;
        core::AtomGroup a = context().groupAtoms(p.group_a);
        core::AtomGroup b = context().groupAtoms(p.group_b);
        if (p.heavy_only) {
            a = core::heavy(a);
            b = core::heavy(b);
        }
        // set(a.ix) == set(b.ix)
        m_same = std::set<std::int32_t>(a.indices().begin(), a.indices().end()) ==
                 std::set<std::int32_t>(b.indices().begin(), b.indices().end());
        if (!m_same) {
            b = without(b, a);
            if (b.empty()) {
                throw core::ValueError("contacts: group " + core::pythonRepr(p.group_b) + " lies entirely inside " +
                                       core::pythonRepr(p.group_a));
            }
        }
        m_a = std::move(a);
        m_b = std::move(b);
        ResidueRows ra = residueRows(m_a);
        ResidueRows rb = residueRows(m_b);
        m_resA = std::move(ra.residues);
        m_rowA = std::move(ra.rows);
        m_resB = std::move(rb.residues);
        m_rowB = std::move(rb.rows);
        m_nB = static_cast<std::int64_t>(m_resB.size());
        m_cutoff = p.cutoff;
        m_atomPairs.assign(frameCount(), 0);
        m_residuePairs.assign(frameCount(), 0);
        m_csr = core::CsrBuilder();
        m_groupA = p.group_a;
        m_groupB = p.group_b;
    }

    void process(std::size_t k, const core::Frame& frame) override {
        const std::vector<float> pa = gather(frame, m_a.indices());
        const std::vector<float> pb = gather(frame, m_b.indices());
        PairList pairs = pairsWithin(pa, pb, m_cutoff, usableBox(context().box(frame)));
        std::vector<std::int64_t>& i = pairs.i;
        std::vector<std::int64_t>& j = pairs.j;
        if (m_same && !i.empty()) {
            // keep a.ix[i] < b.ix[j] and |resindex difference| > 1
            const core::Topology& t = context().topology();
            std::size_t kept = 0;
            for (std::size_t n = 0; n < i.size(); ++n) {
                const std::int32_t ai = m_a[static_cast<std::size_t>(i[n])];
                const std::int32_t bj = m_b[static_cast<std::size_t>(j[n])];
                const std::int64_t ri = t.atomResidue[static_cast<std::size_t>(ai)];
                const std::int64_t rj = t.atomResidue[static_cast<std::size_t>(bj)];
                if (ai < bj && (ri > rj ? ri - rj : rj - ri) > 1) {
                    i[kept] = i[n];
                    j[kept] = j[n];
                    ++kept;
                }
            }
            i.resize(kept);
            j.resize(kept);
        }
        m_atomPairs[k] = static_cast<std::int64_t>(i.size());
        if (i.empty()) {
            m_csr.add({}, {}, {});
            return;
        }
        // np.unique(row_a[i] * n_b + row_b[j], return_counts=True)
        std::vector<std::int64_t> keys(i.size());
        for (std::size_t n = 0; n < i.size(); ++n) {
            keys[n] = m_rowA[static_cast<std::size_t>(i[n])] * m_nB + m_rowB[static_cast<std::size_t>(j[n])];
        }
        std::sort(keys.begin(), keys.end());
        std::vector<std::int64_t> rows, cols;
        std::vector<double> counts;
        for (std::size_t n = 0; n < keys.size();) {
            std::size_t end = n;
            while (end < keys.size() && keys[end] == keys[n]) {
                ++end;
            }
            rows.push_back(keys[n] / m_nB);
            cols.push_back(keys[n] % m_nB);
            counts.push_back(static_cast<double>(end - n));
            n = end;
        }
        m_residuePairs[k] = static_cast<std::int64_t>(rows.size());
        m_csr.add(rows, cols, counts);
    }

    core::AnalysisResult finalize() override {
        const core::Topology& t = context().topology();
        core::AnalysisResult res = newResult({}, {{"labels_a", stringArray(core::residueLabels(t, m_resA))},
                                                  {"labels_b", stringArray(core::residueLabels(t, m_resB))},
                                                  {"group_a", m_groupA},
                                                  {"group_b", m_groupB},
                                                  {"cutoff", m_cutoff}});
        core::NamedArrays arrays = m_csr.arrays("map");
        const auto offsets = arrays[0].second.values<std::int64_t>();
        const auto mi = arrays[1].second.values<std::int32_t>();
        const auto mj = arrays[2].second.values<std::int32_t>();
        const auto mv = arrays[3].second.values<float>();

        // Fraction of native (reference = first analysed frame) residue contacts.
        const std::size_t n = frameCount();
        std::vector<double> q(n, 0.0);
        std::vector<std::int64_t> keys(mi.size());
        for (std::size_t e = 0; e < mi.size(); ++e) {
            keys[e] = static_cast<std::int64_t>(mi[e]) * m_nB + mj[e];
        }
        std::vector<std::int64_t> native;
        if (n > 0) {
            native.assign(keys.begin() + offsets[0], keys.begin() + offsets[1]);
        }
        if (!native.empty()) {
            std::vector<std::int64_t> sortedNative = native;
            std::sort(sortedNative.begin(), sortedNative.end());
            for (std::size_t f = 0; f < n; ++f) {
                for (auto e = static_cast<std::size_t>(offsets[f]); e < static_cast<std::size_t>(offsets[f + 1]); ++e) {
                    if (std::binary_search(sortedNative.begin(), sortedNative.end(), keys[e])) {
                        q[f] += 1.0;
                    }
                }
            }
            for (double& v : q) {
                v /= static_cast<double>(native.size());
            }
        }
        core::NumArray perA = perRowCounts(offsets, mi, mv, m_resA.size());
        core::NumArray perB = perRowCounts(offsets, mj, mv, m_resB.size());
        for (auto& entry : arrays) {
            res.setArray(entry.first, std::move(entry.second));
        }
        res.setArray("atom_pairs", core::NumArray(m_atomPairs));
        res.setArray("residue_pairs", core::NumArray(m_residuePairs));
        res.setArray("q_native", core::NumArray(std::move(q)));
        res.setArray("per_res_a", std::move(perA));
        res.setArray("per_res_b", std::move(perB));
        const std::string& ga = m_groupA;
        const std::string& gb = m_groupB;
        res.addSeries("atom_pairs", "atom_pairs", ga + "–" + gb + " heavy-atom contacts", "count");
        res.addSeries("residue_pairs", "residue_pairs", ga + "–" + gb + " residue contacts", "count");
        res.addSeries("q_native", "q_native", "Fraction of native contacts (Q)", "");
        res.addMatrix("per_residue_a", "per_res_a", "labels_a", ga + " residue contacts with " + gb, "count");
        res.addMatrix("per_residue_b", "per_res_b", "labels_b", gb + " residue contacts with " + ga, "count");
        res.addPairmap("map", "map", "labels_a", "labels_b", ga + "–" + gb + " contact map");
        return res;
    }

private:
    bool m_same = false;
    core::AtomGroup m_a, m_b;
    std::vector<std::int32_t> m_resA, m_resB;
    std::vector<std::int64_t> m_rowA, m_rowB;
    std::int64_t m_nB = 0;
    double m_cutoff = 0.0;
    std::vector<std::int64_t> m_atomPairs, m_residuePairs;
    core::CsrBuilder m_csr;
    std::string m_groupA, m_groupB;
};

class ContactsDescriptor final : public AnalysisDescriptor {
public:
    ContactsDescriptor()
        : AnalysisDescriptor({"contacts",
                              "Pairwise contacts",
                              1,
                              false,
                              2,
                              {{"atom_pairs", "{a}–{b} heavy-atom contacts"},
                               {"residue_pairs", "{a}–{b} residue contacts"},
                               {"q_native", "{a}–{b} fraction of native contacts (Q)"},
                               {"per_residue_a", "{a} residue contacts with {b}"},
                               {"per_residue_b", "{b} residue contacts with {a}"},
                               {"map", "{a}–{b} contact map"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.contacts.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        return detail::sectionParams(project, "contacts");
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"atom_pairs", "residue_pairs", "q_native"}, {"per_residue_a", "per_residue_b"}, {"map"}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::ContactParams& c = project.analyses.contacts;
        return {c.group_a, c.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<ContactsAnalysis>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& contactsAnalysis() {
    static const ContactsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
