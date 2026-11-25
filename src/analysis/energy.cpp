// Non-bonded interaction energy between two groups (mdatk/analysis/energy.py).
//
// E = sum over atom pairs (a in A, b in B, r <= cutoff, not excluded) of
//     332.0522 q_a q_b / (eps r)   (eps r^2 with a distance-dependent dielectric)
//   + A_ab / r^12 - B_ab / r^6
// with Lennard-Jones coefficients from the parm7 topology (use_parm7) or, without them, from Amber
// residue templates and Amber atom types (energy_parameters.cpp), and the magnitude of the net
// force of B on A.
//
// Compiled without floating-point contraction (see CMakeLists.txt): every term is rounded as numpy
// rounds it (the same operations in the same order; integer powers through libm's pow as numpy's
// generic power loop does). The sums over pairs differ from MDATK's in the last bits only because
// the pairs come in another order (pairsWithin sorts them by (i, j); scipy's tree order is not
// reproduced).
//
// Deliberate differences (docs/DIFFERENCES.md): use_parm7 with a topology that is not a parm7 file
// is reported (AN-06); without topology bonds the exclusions follow bonds guessed on the first
// analysed frame, with a warning (AN-07); exclusions follow any bond path (AN-08); without parm7
// parameters, Lennard-Jones parameters (and, without topology charges, partial charges) come from
// Amber force-field residue templates, else Amber atom types, recorded per layer in the metadata
// and cited in the log (AN-25 - AN-27, CIT-04; MDATK: one ad-hoc value per element).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/common.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/detail/numpy_sum.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/io/parm7.hpp"

namespace mdbridge::analysis {

namespace {

core::Array stringArray(const std::vector<std::string>& items) { return core::Array(items.begin(), items.end()); }

// Metadata "lj_table" of a run without parm7 parameters: the name of the table the parameters come
// from, as MDATK 0.1.1 names its generic table (generic_lj.py TABLE_NAME); docs/GENERIC_LJ.md
// records every number of both layers.
constexpr std::string_view kLjTable =
    "Amber force-field residue templates (OpenMM amber14: ff14SB, OL15, OL3, TIP3P with Joung-Cheatham and "
    "Li-Merz ions), then Amber parm10 atom types, Joung-Cheatham and Li-Merz (CM) ions for TIP3P "
    "(docs/GENERIC_LJ.md)";

std::string atomCount(std::int64_t n) { return std::to_string(n) + (n == 1 ? " atom" : " atoms"); }

std::vector<double> widenedPositions(const core::AtomGroup& atoms, const core::Frame& frame) {
    std::vector<double> out;
    out.reserve(atoms.size() * 3);
    for (const std::int32_t atom : atoms.indices()) {
        const auto base = static_cast<std::size_t>(atom) * 3;
        out.push_back(frame.positions.at(base));
        out.push_back(frame.positions.at(base + 1));
        out.push_back(frame.positions.at(base + 2));
    }
    return out;
}

class EnergyAnalysis final : public Analysis {
public:
    EnergyAnalysis(const AnalysisDescriptor& descriptor, AnalysisContext& context, std::string instance)
        : Analysis(descriptor, context, std::move(instance)), m_p(context.project().analyses.energy) {}

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const AnalysisContext& ctx = context();
        m_a = ctx.groupAtoms(m_p.group_a);
        m_b = without(ctx.groupAtoms(m_p.group_b), m_a);
        if (m_b.empty()) {
            throw core::ValueError("energy: group " + core::pythonRepr(core::Value(m_p.group_b)) +
                                   " lies entirely inside " + core::pythonRepr(core::Value(m_p.group_a)));
        }
        const core::Topology& topology = ctx.topology();
        const bool topologyCharges = topology.charges.has_value();
        if (topologyCharges) {
            for (const std::int32_t atom : m_a.indices()) {
                m_qa.push_back((*topology.charges)[static_cast<std::size_t>(atom)]);
            }
            for (const std::int32_t atom : m_b.indices()) {
                m_qb.push_back((*topology.charges)[static_cast<std::size_t>(atom)]);
            }
            m_chargeSource = "topology";
        }
        const std::string& path = ctx.topologyPath();
        if (m_p.use_parm7 && !path.empty() && !io::isParm7(path)) {
            // AN-06: say so (MDATK falls back to generic parameters silently)
            ctx.warn("energy: use_parm7 is set but the topology " + path +
                     " is not a parm7 file; Lennard-Jones parameters come from force-field residue templates and "
                     "Amber atom types");
        }
        if (m_p.use_parm7 && !path.empty() && io::isParm7(path)) {
            try {
                io::LjTable table = io::LjTable::fromParm7(path);
                if (table.typeIndex.size() != topology.atomCount()) {
                    throw core::ValueError("atom count mismatch");
                }
                m_lj = std::move(table);
                m_ljSource = "parm7";
            } catch (const std::exception& e) {
                ctx.warn("energy: could not read LJ parameters from parm7 (" + std::string(e.what()) +
                         "); Lennard-Jones parameters come from force-field residue templates and Amber atom types");
                m_lj.reset();
            }
        }
        const std::size_t na = m_a.size();
        const std::size_t nb = m_b.size();
        m_ljAtoms = {0, 0, 0};
        m_chargeAtoms = {topologyCharges ? static_cast<std::int64_t>(na + nb) : 0, 0, 0};
        m_ljReferences.clear();
        m_defaultElements.clear();
        m_parametersGuessedBonds = false;
        if (m_lj) {
            m_ljAtoms[0] = static_cast<std::int64_t>(na + nb);
        } else {
            // AN-25 - AN-27: templates, then Amber types; template charges without topology charges
            std::vector<std::int32_t> both(m_a.indices());
            both.insert(both.end(), m_b.indices().begin(), m_b.indices().end());
            const core::Frame* frame = ctx.system().currentFrame();
            if (frame == nullptr) {
                throw core::ValueError("energy: no frame to assign parameters on");
            }
            const EnergyParameters params = energyParameters(core::AtomGroup(m_a.sharedTopology(), both), *frame);
            m_ra.assign(params.rstar.begin(), params.rstar.begin() + static_cast<std::ptrdiff_t>(na));
            m_rb.assign(params.rstar.begin() + static_cast<std::ptrdiff_t>(na), params.rstar.end());
            m_ea.assign(params.epsilon.begin(), params.epsilon.begin() + static_cast<std::ptrdiff_t>(na));
            m_eb.assign(params.epsilon.begin() + static_cast<std::ptrdiff_t>(na), params.epsilon.end());
            for (const EnergyParameters::Layer layer : params.layer) {
                ++m_ljAtoms[layer == EnergyParameters::Layer::Template ? 1 : 2];
            }
            m_ljSource = m_ljAtoms[1] == 0 ? "typed" : m_ljAtoms[2] == 0 ? "templates" : "templates+typed";
            m_ljReferences = params.references;
            m_parametersGuessedBonds = params.guessedBonds;
            m_defaultElements = params.defaultElements;
            std::size_t charged = 0;
            if (!topologyCharges) {
                m_qa.assign(params.charge.begin(), params.charge.begin() + static_cast<std::ptrdiff_t>(na));
                m_qb.assign(params.charge.begin() + static_cast<std::ptrdiff_t>(na), params.charge.end());
                charged = static_cast<std::size_t>(
                    std::count(params.templateCharge.begin(), params.templateCharge.end(), true));
                m_chargeAtoms = {0, static_cast<std::int64_t>(charged), static_cast<std::int64_t>(na + nb - charged)};
            }
            if (!topologyCharges && charged == 0) {
                ctx.warn("energy: topology has no partial charges; electrostatics will be zero");
            } else if (!topologyCharges && charged < na + nb) {
                // AN-27: never a partial charge set of a residue; say which atoms have none
                const std::size_t missing = na + nb - charged;
                ctx.warn("energy: the topology has no partial charges; " + std::to_string(charged) + " of " +
                         std::to_string(na + nb) +
                         " atoms take the charges of the force-field residue template they match completely "
                         "(hydrogens included); the other " + atomCount(static_cast<std::int64_t>(missing)) + " (" +
                         std::to_string(params.chargelessResidues) +
                         (params.chargelessResidues == 1 ? " residue" : " residues") +
                         ") have no charge, so their electrostatic terms are zero");
            }
            if (!params.defaultElements.empty()) {
                std::string list;
                for (const std::string& e : params.defaultElements) {
                    list += (list.empty() ? "" : ", ") + core::pythonRepr(core::Value(e));
                }
                ctx.warn("energy: no Lennard-Jones parameters for element(s) " + list + " (" +
                         std::to_string(params.defaultAtoms) + (params.defaultAtoms == 1 ? " atom" : " atoms") +
                         "); carbon's (Amber type C) are used");
            }
            // CIT-04: one line naming what was used and whom to cite
            std::string line = "energy: Lennard-Jones parameters without parm7: ";
            if (m_ljAtoms[1] > 0) {
                std::string ffs;
                for (const std::string& f : params.forceFields) {
                    ffs += (ffs.empty() ? "" : ", ") + f;
                }
                line += atomCount(m_ljAtoms[1]) + " from Amber force-field residue templates (" + ffs + ")";
            }
            if (m_ljAtoms[2] > 0) {
                line += std::string(m_ljAtoms[1] > 0 ? ", " : "") + atomCount(m_ljAtoms[2]) +
                        " from Amber atom types by element";
            }
            line += topologyCharges ? "; partial charges from the topology"
                    : charged > 0   ? "; partial charges of " + atomCount(static_cast<std::int64_t>(charged)) +
                                        " from the templates"
                                    : "; no partial charges";
            line += ". Cite: ";
            for (std::size_t k = 0; k < params.references.size(); ++k) {
                line += (k ? "; " : "") + energyShortCitation(params.references[k]);
            }
            ctx.reporter().log(line);
        }
        if (!topologyCharges) {
            const bool any = m_chargeAtoms[1] > 0;
            if (!any) {
                if (m_lj) {  // parm7 parameters without charges (not written by Amber's tools)
                    ctx.warn("energy: topology has no partial charges; electrostatics will be zero");
                }
                m_qa.assign(na, 0.0);
                m_qb.assign(nb, 0.0);
                m_chargeAtoms = {0, 0, static_cast<std::int64_t>(na + nb)};
            }
            m_chargeSource = !any ? "none" : m_chargeAtoms[2] == 0 ? "templates" : "partial";
        }
        m_excluded.clear();
        m_guessExclusions = m_p.exclude_bonded > 0 && (!topology.bonds || topology.bonds->empty());
        if (m_guessExclusions) {
            // AN-07: MDATK excludes nothing, silently
            // CIT-03: the guess uses Cordero et al.'s covalent radii, cited in the warning
            ctx.warn("energy: the topology has no bonds; bonded pairs are excluded along bonds guessed from "
                     "distances on the first analysed frame (covalent radii: Cordero et al., Dalton Trans. "
                     "2008, 2832-2838, doi:10.1039/b801115j)");
        } else {
            m_excluded = energyBondedExclusions(m_a, m_b, m_p.exclude_bonded);
        }
        const std::size_t n = frameCount();
        m_elec.assign(n, 0.0);
        m_vdw.assign(n, 0.0);
        m_force.assign(n, 0.0);
        m_pairCount.assign(n, 0);
    }

    // AN-07: without topology bonds, the exclusions follow the bonds guessed on the first analysed
    // frame (guessBonds over the whole system, so paths through other atoms count too)
    // CIT-06: the references this run used: those of its parameters (lj_references), then the
    // covalent radii (Cordero et al. 2008) when bonds were guessed (exclusions, AN-07; parameters)
    std::vector<std::string> citationKeys() const {
        std::vector<std::string> keys = m_ljReferences;
        if (m_guessExclusions || m_parametersGuessedBonds) {
            keys.emplace_back("cordero2008");
        }
        return keys;
    }

    void prepare(const core::Frame& first) override {
        if (!m_guessExclusions) {
            return;
        }
        const std::vector<std::int64_t> guessed = guessBonds(context().system().atoms(), first);
        std::vector<core::Bond> bonds;
        bonds.reserve(guessed.size() / 2);
        for (std::size_t p = 0; p + 1 < guessed.size(); p += 2) {
            bonds.push_back({static_cast<std::int32_t>(guessed[p]), static_cast<std::int32_t>(guessed[p + 1])});
        }
        m_excluded = energyBondedExclusions(m_a, m_b, m_p.exclude_bonded, &bonds);
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        const std::optional<Box> box = usableBox(context().box(frame));
        const std::vector<double> pa = widenedPositions(m_a, frame);
        const std::vector<double> pb = widenedPositions(m_b, frame);
        PairList pairs = pairsWithin(pa, pb, m_p.cutoff, box);
        if (!m_excluded.empty() && !pairs.empty()) {
            const auto nb = static_cast<std::int64_t>(m_b.size());
            std::size_t kept = 0;
            for (std::size_t p = 0; p < pairs.size(); ++p) {
                if (!std::binary_search(m_excluded.begin(), m_excluded.end(), pairs.i[p] * nb + pairs.j[p])) {
                    pairs.i[kept] = pairs.i[p];
                    pairs.j[kept] = pairs.j[p];
                    pairs.d[kept] = pairs.d[p];
                    ++kept;
                }
            }
            pairs.i.resize(kept);
            pairs.j.resize(kept);
            pairs.d.resize(kept);
        }
        const std::size_t n = pairs.size();
        m_pairCount[k] = static_cast<std::int64_t>(n);
        if (n == 0) {
            return;
        }
        std::vector<std::int64_t> atomI, atomJ;
        if (m_lj) {
            atomI.reserve(n);
            atomJ.reserve(n);
            for (std::size_t p = 0; p < n; ++p) {
                atomI.push_back(m_a[static_cast<std::size_t>(pairs.i[p])]);
                atomJ.push_back(m_b[static_cast<std::size_t>(pairs.j[p])]);
            }
        }
        std::vector<double> acoef, bcoef;
        if (m_lj) {
            std::tie(acoef, bcoef) = m_lj->pairCoefficients(atomI, atomJ);
        }
        std::vector<double> eElec(n), eVdw(n), fmag(n);
        for (std::size_t p = 0; p < n; ++p) {
            const auto i = static_cast<std::size_t>(pairs.i[p]);
            const auto j = static_cast<std::size_t>(pairs.j[p]);
            const double d = std::max(pairs.d[p], 0.5);  // np.maximum(d, 0.5)
            const double qq = kEnergyCoulomb * m_qa[i] * m_qb[j] / m_p.dielectric;
            double deElec = 0.0;
            if (m_p.distance_dependent) {
                eElec[p] = qq / (d * d);
                deElec = -2.0 * qq / std::pow(d, 3.0);
            } else {
                eElec[p] = qq / d;
                deElec = -qq / (d * d);
            }
            double a = 0.0;
            double b = 0.0;
            if (m_lj) {
                a = acoef[p];
                b = bcoef[p];
            } else {
                const double rmin = m_ra[i] + m_rb[j];
                const double eps = std::sqrt(m_ea[i] * m_eb[j]);
                const double r6 = std::pow(rmin, 6.0);
                a = eps * r6 * r6;
                b = 2.0 * eps * r6;
            }
            const double inv6 = 1.0 / std::pow(d, 6.0);
            eVdw[p] = a * inv6 * inv6 - b * inv6;
            const double deVdw = (-12.0 * a * inv6 * inv6 + 6.0 * b * inv6) / d;
            // the scalar force at the clamped distance along the unit vector of the real one:
            // rvec * (F / r) (MDATK divides by the clamped d, which shrinks the force by r / 0.5
            // for r < 0.5; docs/DIFFERENCES.md AN-05); no direction for r = 0
            const double r = pairs.d[p];
            fmag[p] = r > 0 ? -(deElec + deVdw) / r : 0.0;
        }
        m_elec[k] = core::detail::numpySum(std::span<const double>(eElec));
        m_vdw[k] = core::detail::numpySum(std::span<const double>(eVdw));
        // rvec = minimum_image_vectors(pa[i] - pb[j], box); (rvec * fmag[:, None]).sum(axis=0)
        std::vector<double> rvec(3 * n);
        for (std::size_t p = 0; p < n; ++p) {
            const auto i = static_cast<std::size_t>(pairs.i[p]);
            const auto j = static_cast<std::size_t>(pairs.j[p]);
            for (std::size_t c = 0; c < 3; ++c) {
                rvec[3 * p + c] = pa[3 * i + c] - pb[3 * j + c];
            }
        }
        if (box) {
            rvec = minimumImageVectors(rvec, box);  // the image the pair search measured
        }
        double total[3] = {0.0, 0.0, 0.0};
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t c = 0; c < 3; ++c) {
                total[c] += rvec[3 * p + c] * fmag[p];
            }
        }
        // np.linalg.norm of a 3-vector: sqrt(x . x)
        m_force[k] = std::sqrt((total[0] * total[0] + total[1] * total[1]) + total[2] * total[2]);
    }

    core::AnalysisResult finalize() override {
        core::AnalysisResult res = newResult({}, {{"lj_source", m_ljSource},
                                                  {"charge_source", m_chargeSource},
                                                  {"group_a", m_p.group_a},
                                                  {"group_b", m_p.group_b},
                                                  {"n_excluded_pairs", static_cast<std::int64_t>(m_excluded.size())},
                                                  {"lj_atoms", core::Table{{"parm7", m_ljAtoms[0]},
                                                                           {"template", m_ljAtoms[1]},
                                                                           {"typed", m_ljAtoms[2]}}},
                                                  {"charge_atoms", core::Table{{"topology", m_chargeAtoms[0]},
                                                                               {"template", m_chargeAtoms[1]},
                                                                               {"none", m_chargeAtoms[2]}}},
                                                  {"lj_references", stringArray(m_ljReferences)},
                                                  {"lj_default_elements", stringArray(m_defaultElements)},
                                                  {"citation_keys", stringArray(citationKeys())}});
        if (m_ljSource != "parm7") {
            // MDATK 0.1.1's key, after n_excluded_pairs as MDATK has it (a run without parm7 parameters)
            for (auto it = res.meta.begin(); it != res.meta.end(); ++it) {
                if (it->first == "n_excluded_pairs") {
                    res.meta.insert(it + 1, {"lj_table", core::Value(std::string(kLjTable))});
                    break;
                }
            }
        }
        std::vector<double> total(m_elec.size());
        for (std::size_t k = 0; k < total.size(); ++k) {
            total[k] = m_elec[k] + m_vdw[k];
        }
        res.arrays.emplace_back("total", core::NumArray(std::move(total)));
        res.arrays.emplace_back("elec", core::NumArray(m_elec));
        res.arrays.emplace_back("vdw", core::NumArray(m_vdw));
        res.arrays.emplace_back("force", core::NumArray(m_force));
        res.arrays.emplace_back("n_pairs", core::NumArray(m_pairCount));
        const std::string& ga = m_p.group_a;
        const std::string& gb = m_p.group_b;
        const core::Value target = m_p.target_value ? core::Value(*m_p.target_value) : core::Value();
        res.addSeries("total", "total", ga + "–" + gb + " interaction energy", "kcal/mol", {{"target", target}});
        res.addSeries("elec", "elec", ga + "–" + gb + " electrostatic energy", "kcal/mol");
        res.addSeries("vdw", "vdw", ga + "–" + gb + " van der Waals energy", "kcal/mol");
        res.addSeries("net_force", "force", "Net force of " + gb + " on " + ga, "kcal/mol/Å");
        return res;
    }

private:
    core::EnergyParams m_p;
    core::AtomGroup m_a;
    core::AtomGroup m_b;
    std::vector<double> m_qa, m_qb;
    std::string m_chargeSource;
    std::string m_ljSource;
    std::optional<io::LjTable> m_lj;
    std::vector<double> m_ra, m_ea, m_rb, m_eb;
    std::array<std::int64_t, 3> m_ljAtoms{};      // parm7, template, typed
    std::array<std::int64_t, 3> m_chargeAtoms{};  // topology, template, none
    std::vector<std::string> m_ljReferences;
    std::vector<std::string> m_defaultElements;
    bool m_parametersGuessedBonds = false;
    std::vector<std::int64_t> m_excluded;
    bool m_guessExclusions = false;
    std::vector<double> m_elec, m_vdw, m_force;
    std::vector<std::int64_t> m_pairCount;
};

class EnergyDescriptor final : public AnalysisDescriptor {
public:
    EnergyDescriptor()
        : AnalysisDescriptor({"energy",
                              "Intermolecular energy",
                              2,  // v2: exclusions along any bond path, guessed bonds (AN-07, AN-08);
                                  // float64 parm7 charges (RD-P12); periodic geometry and force
                                  // direction (AN-01-AN-05); per-atom element fallback (SEL-13);
                                  // template / typed parameters and template charges (AN-25-AN-27)
                              false,
                              2,
                              {{"total", "{a}–{b} interaction energy"},
                               {"elec", "{a}–{b} electrostatic energy"},
                               {"vdw", "{a}–{b} van der Waals energy"},
                               {"net_force", "Net force of {b} on {a}"}}}) {}

    bool enabled(const core::Project& project) const override { return project.analyses.energy.enabled; }

    core::Table paramsFor(const core::Project& project) const override {
        core::Table params = detail::sectionParams(project, "energy");
        // target_value is display only
        std::erase_if(params, [](const auto& entry) { return entry.first == "target_value"; });
        return params;
    }

    ExpectedOutputs expectedOutputs(const core::Project&) const override {
        return {{"total", "elec", "vdw", "net_force"}, {}, {}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        const core::EnergyParams& e = project.analyses.energy;
        return {e.group_a, e.group_b};
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<EnergyAnalysis>(*this, context, std::move(instance));
    }
};

}  // namespace

std::vector<std::int64_t> energyBondedExclusions(const core::AtomGroup& a, const core::AtomGroup& b,
                                                 std::int64_t depth, const std::vector<core::Bond>* bonds) {
    std::vector<std::int64_t> keys;
    if (depth <= 0 || a.sharedTopology() == nullptr) {
        return keys;
    }
    const core::Topology& topology = a.topology();
    if (!bonds) {
        if (!topology.bonds) {
            return keys;
        }
        bonds = &*topology.bonds;
    }
    if (bonds->empty()) {  // no bond information: nothing excluded
        return keys;
    }
    const std::size_t nAtoms = topology.atomCount();
    const auto ix = [](std::int32_t atom) { return static_cast<std::size_t>(atom); };
    std::vector<std::vector<std::int32_t>> adjacency(nAtoms);
    for (const core::Bond& bond : *bonds) {
        adjacency[ix(bond[0])].push_back(bond[1]);
        adjacency[ix(bond[1])].push_back(bond[0]);
    }
    // position in b (the last one for a repeated atom, as MDATK's dicts keep it)
    std::vector<std::int64_t> bLocal(nAtoms, -1);
    for (std::size_t n = 0; n < b.size(); ++n) {
        bLocal[ix(b[n])] = static_cast<std::int64_t>(n);
    }
    const auto nb = static_cast<std::int64_t>(b.size());
    // AN-08: from every A atom, breadth-first through any atoms up to `depth` bonds; every B atom
    // reached is excluded (MDATK: only paths starting with a direct A-B bond's A side, through A)
    std::vector<std::int64_t> seen(nAtoms, -1);
    std::vector<std::int32_t> visited;
    std::vector<std::int32_t> queue;
    for (std::size_t n = 0; n < a.size(); ++n) {
        const std::int32_t start = a[n];
        for (const std::int32_t atom : visited) {
            seen[ix(atom)] = -1;
        }
        visited.assign(1, start);
        seen[ix(start)] = 0;
        queue.assign(1, start);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const std::int32_t cur = queue[head];
            if (seen[ix(cur)] >= depth) {
                continue;
            }
            for (const std::int32_t next : adjacency[ix(cur)]) {
                if (seen[ix(next)] >= 0) {
                    continue;
                }
                seen[ix(next)] = seen[ix(cur)] + 1;
                visited.push_back(next);
                queue.push_back(next);
                if (bLocal[ix(next)] >= 0) {
                    keys.push_back(static_cast<std::int64_t>(n) * nb + bLocal[ix(next)]);
                }
            }
        }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

const AnalysisDescriptor& energyAnalysis() {
    static const EnergyDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
