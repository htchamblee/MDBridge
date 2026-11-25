// Root mean square deviation against a reference frame (mdatk/analysis/rmsd.py).
//
// Per target: without superposition the plain RMSD to the reference positions; with
// superposition and no fit selection MDAnalysis' rms.rmsd (superposedRmsd, bit for bit); with a
// fit selection the fit atoms are superposed with kabsch and the target atoms moved with that
// rotation. setup() moves the trajectory to the reference frame (MDATK seeks it there), so the
// analyses set up after this one see the reference frame as the current frame.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/analysis/rmsd_qcp.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::analysis {

std::string rmsdOutputKey(std::string_view label, std::size_t k) {
    // re.sub(r"[^A-Za-z0-9_.+\-]+", "_", label).strip("_")
    const auto keep = [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
               c == '+' || c == '-';
    };
    std::string key;
    bool inRun = false;
    for (const char c : label) {
        if (keep(c)) {
            key += c;
            inRun = false;
        } else if (!inRun) {
            key += '_';
            inRun = true;
        }
    }
    const std::size_t first = key.find_first_not_of('_');
    if (first == std::string::npos) {
        return "target" + std::to_string(k);
    }
    return key.substr(first, key.find_last_not_of('_') - first + 1);
}

namespace {

// One target of a run (MDATK's per-target dict).
struct RmsdTarget {
    std::string label;  // t.label or "target<k>"
    std::string out;  // _out_key(t.label, k)
    core::AtomGroup mob;
    std::optional<core::AtomGroup> fit;  // a non-blank fit selection
    std::optional<std::vector<double>> w;  // masses of mob when mass weighted
    std::optional<std::vector<double>> wf;  // masses of fit (w without one) when mass weighted
    std::vector<double> ref;  // mob positions in the reference frame (float64, n x 3)
    std::vector<Vec3> refFit;  // fit positions in the reference frame
    bool translationOnly = false;  // AN-18: fewer than 3 or collinear fit atoms
};

// AN-18: whether the fit atoms (reference positions, weights or none) are fewer than 3 or
// collinear: q <= 1e-10 p^2 for the weighted inner-product matrix of the centred set (p its trace,
// q the sum of its principal 2 x 2 minors), the rule of the superposition RMSD (rmsd_qcp.cpp).
bool collinearFit(const std::vector<Vec3>& x, const std::optional<std::vector<double>>& w) {
    const std::size_t n = x.size();
    if (n < 3) {
        return true;
    }
    double total = 0.0;
    Vec3 centre{};
    for (std::size_t i = 0; i < n; ++i) {
        const double wi = w ? (*w)[i] : 1.0;
        total += wi;
        for (std::size_t c = 0; c < 3; ++c) {
            centre[c] += wi * x[i][c];
        }
    }
    if (!(total > 0.0)) {
        return false;  // left to kabsch (MDATK's behaviour)
    }
    for (double& c : centre) {
        c /= total;
    }
    double m[3][3] = {};
    for (std::size_t i = 0; i < n; ++i) {
        const double wi = w ? (*w)[i] : 1.0;
        const Vec3 d{x[i][0] - centre[0], x[i][1] - centre[1], x[i][2] - centre[2]};
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                m[r][c] += wi * (d[r] * d[c]);
            }
        }
    }
    const double p = m[0][0] + m[1][1] + m[2][2];
    const double q = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) + (m[0][0] * m[2][2] - m[0][2] * m[2][0]) +
                     (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
    return q <= 1e-10 * p * p;
}

// atoms.positions.astype(np.float64) of `frame`, x, y, z interleaved
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

std::vector<Vec3> vectorsOf(const core::AtomGroup& atoms, const core::Frame& frame) {
    std::vector<Vec3> out;
    out.reserve(atoms.size());
    for (const std::int32_t atom : atoms.indices()) {
        const auto base = 3 * static_cast<std::size_t>(atom);
        out.push_back({static_cast<double>(frame.positions.at(base)), static_cast<double>(frame.positions.at(base + 1)),
                       static_cast<double>(frame.positions.at(base + 2))});
    }
    return out;
}

// atoms.masses.astype(np.float64) (MDAnalysis always has masses; the loaders fill them)
std::vector<double> massesOf(const core::AtomGroup& atoms) {
    const auto& masses = atoms.topology().masses;
    if (!masses) {
        throw std::logic_error("RMSD: the topology has no masses");
    }
    std::vector<double> out;
    out.reserve(atoms.size());
    for (const std::int32_t atom : atoms.indices()) {
        out.push_back((*masses)[static_cast<std::size_t>(atom)]);
    }
    return out;
}

class RmsdRun final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::RMSDParams& p = project().analyses.rmsd;
        if (p.targets.empty()) {
            throw core::ValueError("RMSD: no targets defined");
        }
        const auto n = static_cast<std::int64_t>(frameCount());
        m_refPos = std::min<std::int64_t>(std::max<std::int64_t>(0, p.reference_frame), n - 1);
        if (n > 0 && p.reference_frame >= n) {  // MDATK clamps silently (PRJ-02)
            context().warn("RMSD: reference frame " + std::to_string(p.reference_frame) + " is beyond the " +
                           std::to_string(n) + " analysed frames; using the last");
        } else if (n > 0 && p.reference_frame < 0) {
            context().warn("RMSD: reference frame " + std::to_string(p.reference_frame) +
                           " is negative; using the first");
        }
        m_superposition = p.superposition;
        m_targets.clear();
        for (std::size_t k = 0; k < p.targets.size(); ++k) {
            const core::RMSDTarget& t = p.targets[k];
            RmsdTarget target;
            target.mob = context().select(t.selection, true);
            if (target.mob.empty()) {
                context().warn("RMSD target " + core::pythonStringRepr(t.label) + ": " +
                               core::pythonStringRepr(t.selection) + " matched no atoms; its series will be empty (NaN)");
            }
            if (!core::pythonStrip(t.fit_selection).empty()) {
                target.fit = context().select(t.fit_selection);
            }
            if (p.mass_weighted) {
                target.w = massesOf(target.mob);
                target.wf = target.fit ? massesOf(*target.fit) : target.w;
            }
            target.label = t.label.empty() ? "target" + std::to_string(k) : t.label;
            target.out = rmsdOutputKey(t.label, k);
            m_targets.push_back(std::move(target));
        }
        // self.ctx.u.trajectory[int(frame_indices[self.ref_pos])]: the reference frame stays the
        // current frame for the setups that follow
        if (m_refPos < 0) {
            throw std::out_of_range("index -1 is out of bounds for axis 0 with size 0");
        }
        core::Trajectory* trajectory = context().system().trajectory.get();
        if (!trajectory) {
            throw std::logic_error("RMSD: the system has no trajectory");
        }
        const core::Frame& reference = trajectory->seek(frames[static_cast<std::size_t>(m_refPos)]);
        for (std::size_t k = 0; k < m_targets.size(); ++k) {
            RmsdTarget& t = m_targets[k];
            t.ref = positionsOf(t.mob, reference);
            if (t.fit) {
                t.refFit = vectorsOf(*t.fit, reference);
                t.translationOnly = m_superposition && collinearFit(t.refFit, t.wf);
                if (t.translationOnly) {
                    const std::size_t atoms = t.fit->size();
                    context().warn("RMSD target " + core::pythonStringRepr(p.targets[k].label) + ": fit selection has " +
                                   std::to_string(atoms) + (atoms < 3 ? (atoms == 1 ? " atom" : " atoms")
                                                                      : " collinear atoms") +
                                   "; superposition is translation only");
                }
            }
        }
        m_values.assign(m_targets.size() * frameCount(), std::numeric_limits<double>::quiet_NaN());
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        for (std::size_t n = 0; n < m_targets.size(); ++n) {
            const RmsdTarget& t = m_targets[n];
            if (t.mob.empty()) {
                continue;  // the row stays NaN
            }
            double& value = m_values[n * frameCount() + k];
            std::vector<double> x = positionsOf(t.mob, frame);
            if (m_superposition && !t.fit) {
                // rms.rmsd(x, ref, weights=w, center=True, superposition=True)
                value = t.w ? superposedRmsd(x, t.ref, std::span<const double>(*t.w)) : superposedRmsd(x, t.ref);
                continue;
            }
            const std::size_t atoms = t.mob.size();
            if (!m_superposition) {
                for (std::size_t i = 0; i < x.size(); ++i) {
                    x[i] = x[i] - t.ref[i];  // diff = x - ref
                }
            } else {
                const std::vector<Vec3> xf = vectorsOf(*t.fit, frame);
                const KabschResult fit =
                    t.wf ? kabsch(xf, t.refFit, std::span<const double>(*t.wf)) : kabsch(xf, t.refFit);
                // AN-18: translation only for fewer than 3 or collinear fit atoms (MDATK: an
                // arbitrary rotation about the line for collinear ones)
                static constexpr Mat3 kIdentity{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
                const Mat3& r = t.translationOnly ? kIdentity : fit.rotation;
                // diff = (x - mc) @ r.T + rc - ref
                for (std::size_t i = 0; i < atoms; ++i) {
                    const double a0 = x[3 * i] - fit.mobileCenter[0];
                    const double a1 = x[3 * i + 1] - fit.mobileCenter[1];
                    const double a2 = x[3 * i + 2] - fit.mobileCenter[2];
                    for (std::size_t c = 0; c < 3; ++c) {
                        const double moved = ((a0 * r[c][0] + a1 * r[c][1]) + a2 * r[c][2]) + fit.referenceCenter[c];
                        x[3 * i + c] = moved - t.ref[3 * i + c];
                    }
                }
            }
            // sq = (diff * diff).sum(axis=1)
            std::vector<double> sq(atoms);
            for (std::size_t i = 0; i < atoms; ++i) {
                const double* d = &x[3 * i];
                sq[i] = core::detail::numpySum<double>([d](std::size_t c) { return d[c] * d[c]; }, 3);
            }
            if (t.w) {
                // np.average(sq, weights=w): sum(sq * w) / sum(w)
                const std::vector<double>& w = *t.w;
                const double scale = core::detail::numpySum(w);
                if (scale == 0.0) {
                    throw ZeroDivisionError("Weights sum to zero, can't be normalized");
                }
                const double total =
                    core::detail::numpySum<double>([&sq, &w](std::size_t i) { return sq[i] * w[i]; }, atoms);
                value = std::sqrt(total / scale);
            } else {
                value = std::sqrt(core::detail::numpyMean(sq));  // np.sqrt(sq.mean())
            }
        }
    }

    core::AnalysisResult finalize() override {
        core::AnalysisResult res = newResult({}, core::Table{{"reference_frame_index", m_refPos}});
        const std::vector<std::string> groups = parseInstance(instance()).groups;
        const std::string suffix = groups.empty() ? std::string() : " (" + groups[0] + ")";
        const std::size_t frames = frameCount();
        for (std::size_t n = 0; n < m_targets.size(); ++n) {
            const std::string arr = "rmsd_" + std::to_string(n);
            const auto first = m_values.begin() + static_cast<std::ptrdiff_t>(n * frames);
            res.setArray(arr, core::NumArray(std::vector<double>(first, first + static_cast<std::ptrdiff_t>(frames))));
            res.addSeries(m_targets[n].out, arr, "RMSD " + m_targets[n].label + suffix, "\xC3\x85");  // Å
        }
        return res;
    }

private:
    std::int64_t m_refPos = 0;
    bool m_superposition = true;
    std::vector<RmsdTarget> m_targets;
    std::vector<double> m_values;  // targets x frames, NaN where not computed
};

class RmsdDescriptor final : public AnalysisDescriptor {
public:
    // version 2: exact values for degenerate sets and near 0 (AN-16, AN-17); translation only for
    // fewer than 3 or collinear fit atoms (AN-18)
    RmsdDescriptor() : AnalysisDescriptor({"rmsd", "RMSD", 2, false, 1, {}}) {}

    void withGroups(core::Project& project, const std::vector<std::string>& groups) const override {
        const std::string& g = groups.at(0);
        project.analyses.rmsd.targets = {{"Heavy", "group:" + g + " and not name H*", ""},
                                         {"Backbone", "group:" + g + " and name C4'", ""}};
    }

    std::string outputLabel(const core::Project& project, std::string_view out) const override {
        const auto& targets = project.analyses.rmsd.targets;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            if (rmsdOutputKey(targets[k].label, k) == out) {
                return "RMSD " + targets[k].label;
            }
        }
        return "RMSD " + std::string(out);
    }

    bool enabled(const core::Project& project) const override {
        return project.analyses.rmsd.enabled && !project.analyses.rmsd.targets.empty();
    }

    core::Table paramsFor(const core::Project& project) const override {
        const core::Table section = detail::sectionParams(project, "rmsd");
        core::Table out;
        for (const char* key : {"targets", "reference_frame", "superposition", "mass_weighted"}) {
            for (const auto& [k, v] : section) {
                if (k == key) {
                    out.emplace_back(k, v);
                }
            }
        }
        return out;
    }

    ExpectedOutputs expectedOutputs(const core::Project& project) const override {
        ExpectedOutputs out;
        const auto& targets = project.analyses.rmsd.targets;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            out.series.push_back(rmsdOutputKey(targets[k].label, k));
        }
        return out;
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out;
        for (const core::RMSDTarget& t : project.analyses.rmsd.targets) {
            out.push_back(t.selection);
            out.push_back(t.fit_selection);
        }
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<RmsdRun>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& rmsdAnalysis() {
    static const RmsdDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
