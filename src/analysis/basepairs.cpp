// Base pairs of the structural-view residues (mdatk/analysis/basepairs.py; internal, feeds the
// secondary-structure panels).
//
// Rows are exactly the residues of the trace analysis (structureResidues, same order). Pairs are
// detected with the builtin detector (chem::PairDetector: canonical Watson-Crick and optional G-U
// wobble pairs, each base in at most one pair) on every stride-th analysed frame and held in
// between. Without use_pbc distances ignore periodic images, as in MDATK. With use_pbc and a
// periodic cell every donor-acceptor and C1'-C1' distance is a minimum-image distance (MDATK found
// no pairs in a duplex wrapped across a face; docs/DIFFERENCES.md AN-11). Per-distance images
// rather than whole residues: the distances read are at most 12.5 A, while residues' reference
// atoms across a pair (P-P about 18 A) can be more than half a small cell apart.
//
// Outputs: pair map "map" (value 1 Watson-Crick, 2 wobble, 3 other; each pair once with i < j),
// series "count" (pairs per frame), array "resindex" and meta residue_labels, engine, stride,
// value_labels.
//
// Deliberate difference: MDBridge does not bundle Barnaba (GPL). With engine = "barnaba" the
// analysis warns and runs the builtin detector instead (meta engine "builtin", no lw_first_frame,
// noncanonical has no effect); params keep the project's values, so the cache key is MDATK's.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/selections.hpp"

namespace mdbridge::analysis {

namespace {

constexpr const char* kBarnabaFallback = "basepairs: Barnaba is not available in MDBridge; using the builtin detector";

class BasePairs final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::BasePairParams& bp = project().analyses.basepairs;
        StructureResidues structure = structureResidues(context());
        if (structure.residues.empty()) {
            throw core::SelectionError("base pairs: no residues selected for the structural view");
        }
        m_residues = std::move(structure.residues);
        m_stride = bp.stride < 1 ? 1 : bp.stride;
        if (bp.engine == "barnaba") {
            context().warn(kBarnabaFallback);
        }
        const core::Topology& t = context().topology();
        const core::AtomGroup all = context().system().atoms();
        std::vector<std::string> resnames;
        std::vector<chem::AtomIndex> atomIndex;
        m_atoms.clear();
        for (const std::int32_t r : m_residues) {
            resnames.push_back(t.resname(r));
            const auto span = t.residueAtoms(r);
            // heavy(res.atoms): elements decided on the residue's own atoms
            const core::AtomGroup residueAtoms(all.sharedTopology(), std::vector<std::int32_t>(span.begin(), span.end()));
            const core::AtomGroup heavyAtoms = core::heavy(residueAtoms);
            chem::AtomIndex names;
            for (const std::int32_t a : heavyAtoms.indices()) {
                std::string name = t.name(a);
                for (char& c : name) {
                    if (c == '*') {
                        c = '\'';  // old-style primes
                    }
                }
                // first alternate location wins
                if (names.emplace(std::move(name), static_cast<std::int64_t>(m_atoms.size())).second) {
                    m_atoms.push_back(a);
                }
            }
            atomIndex.push_back(std::move(names));
        }
        m_detector = std::make_unique<chem::PairDetector>(resnames, atomIndex, bp.cutoff, bp.wobble);
        m_framePairs.assign(frameCount(), {});
        m_source.assign(frameCount(), 0);
    }

    bool concurrentProcess() const override { return true; }

    // Frames between two strides hold the pairs of the last detected frame (k - k % stride).
    void process(std::size_t k, const core::Frame& frame) override {
        const auto offset = static_cast<std::size_t>(static_cast<std::int64_t>(k) % m_stride);
        if (offset != 0) {
            m_source[k] = k - offset;
            return;
        }
        m_source[k] = k;
        // self.atoms.positions.astype(np.float64)
        std::vector<double> xyz;
        xyz.reserve(m_atoms.size() * 3);
        for (const std::int32_t a : m_atoms) {
            const auto at = static_cast<std::size_t>(a) * 3;
            for (std::size_t d = 0; d < 3; ++d) {
                xyz.push_back(static_cast<double>(frame.positions[at + d]));
            }
        }
        // with use_pbc and a periodic cell, every distance the detector reads is a minimum-image
        // distance, so pairs across a face of a wrapped trajectory are found (AN-11)
        if (const std::optional<PeriodicCell> cell = PeriodicCell::from(context().box(frame))) {
            const chem::PairDetector::DifferenceImage image = [&cell](const std::array<double, 3>& v) {
                return cell->minimumImage(v);
            };
            m_framePairs[k] = m_detector->detect(xyz, &image);
            return;
        }
        m_framePairs[k] = m_detector->detect(xyz);
    }

    core::AnalysisResult finalize() override {
        const std::size_t n = frameCount();
        core::CsrBuilder csr;
        std::vector<std::int64_t> count(n, 0);
        for (std::size_t k = 0; k < n; ++k) {
            const std::vector<chem::BasePair>& pairs = m_framePairs[m_source[k]];
            std::vector<std::int64_t> i, j;
            std::vector<double> v;
            for (const chem::BasePair& p : pairs) {
                i.push_back(p.i);
                j.push_back(p.j);
                v.push_back(static_cast<double>(p.kind));
            }
            csr.add(i, j, v);
            count[k] = static_cast<std::int64_t>(pairs.size());
        }
        const core::Topology& t = context().topology();
        const std::vector<std::string> labels = core::residueLabels(t, m_residues);
        core::Table valueLabels;
        for (const int kind : {chem::kWatsonCrick, chem::kWobble, chem::kOther}) {
            valueLabels.emplace_back(std::to_string(kind), std::string(chem::kindLabel(kind)));
        }
        core::AnalysisResult res = newResult({}, {{"residue_labels", core::Array(labels.begin(), labels.end())},
                                                  {"engine", std::string("builtin")},
                                                  {"stride", m_stride},
                                                  {"value_labels", std::move(valueLabels)}});
        for (auto& [name, array] : csr.arrays("map")) {
            res.arrays.emplace_back(name, std::move(array));
        }
        res.arrays.emplace_back("count", core::NumArray(std::move(count)));
        res.arrays.emplace_back("resindex", core::NumArray(std::vector<std::int64_t>(m_residues.begin(), m_residues.end())));
        res.addSeries("count", "count", "Base pairs", "count");
        res.addPairmap("map", "map", "residue_labels", "residue_labels", "Base pairs");
        return res;
    }

private:
    std::vector<std::int32_t> m_residues;
    std::int64_t m_stride = 1;
    std::vector<std::int32_t> m_atoms;  // the rows of the coordinates given to the detector
    std::unique_ptr<chem::PairDetector> m_detector;
    std::vector<std::vector<chem::BasePair>> m_framePairs;
    std::vector<std::size_t> m_source;  // the detection frame each frame copies
};

class BasePairsDescriptor final : public AnalysisDescriptor {
public:
    BasePairsDescriptor()
        : AnalysisDescriptor(
              {"basepairs", "Base pairs", 2, true, 0, {{"map", "Base pairs"}, {"count", "Base pairs"}}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::BasePairParams& bp = project.analyses.basepairs;
        const std::vector<std::string> names = structureGroups(project);
        return {{"groups", core::Array(names.begin(), names.end())},
                {"ions", project.analyses.mg.mg_selection},
                {"engine", bp.engine},
                {"stride", bp.stride},
                {"cutoff", bp.cutoff},
                {"wobble", bp.wobble},
                {"noncanonical", bp.noncanonical}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override { return {{"count"}, {}, {"map"}}; }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<BasePairs>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& basepairsAnalysis() {
    static const BasePairsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
