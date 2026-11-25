// Internal analyses feeding the structural view (mdatk/analysis/structure.py): backbone trace,
// proximity to a residue set, atomic coordinates; and the helpers they share (structure groups and
// residues, group membership, trace points, the bonds of the atoms).
//
// Deliberate difference: when the topology has no bonds among the atoms, the atoms analysis
// guesses them from distances like MDATK, but its `bonds` rows come out sorted by (i, j); MDATK's
// come in the order of scipy's kd-tree traversal (pairs_within without a box), which depends on
// how the tree was built (see spatial.hpp pairsWithin). The set of rows is the same.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "builtin.hpp"
#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/spatial.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::analysis {

std::vector<std::string> structureGroups(const core::Project& project) {
    const std::vector<std::string> all = core::groupNames(project);
    const std::vector<std::string>& shown = project.visual.structure_groups;
    const std::vector<std::string>& names = shown.empty() ? all : shown;
    const std::set<std::string> known(all.begin(), all.end());
    std::vector<std::string> out;
    for (const std::string& n : names) {
        if (known.contains(n)) {
            out.push_back(n);
        }
    }
    return out;
}

std::string structureSelectionText(const core::Project& project) {
    const std::vector<std::string> names = structureGroups(project);
    if (names.empty()) {
        throw core::SelectionError("no residue groups available for the structural view");
    }
    std::string out;
    for (const std::string& n : names) {
        if (!out.empty()) {
            out += " or ";
        }
        out += "group:" + n;
    }
    return out;
}

StructureResidues structureResidues(const AnalysisContext& context) {
    const core::Project& project = context.project();
    core::AtomGroup atoms = context.select(structureSelectionText(project));
    const core::AtomGroup ions = context.select(project.analyses.mg.mg_selection, true);
    if (!ions.empty()) {
        const std::set<std::int32_t> drop(ions.indices().begin(), ions.indices().end());
        std::vector<bool> keep;
        keep.reserve(atoms.size());
        for (const std::int32_t a : atoms.indices()) {
            keep.push_back(!drop.contains(a));
        }
        atoms = atoms.masked(keep);
    }
    std::vector<std::int32_t> residues = atoms.residues();
    return {std::move(residues), std::move(atoms)};
}

std::string structureGuessElement(std::string_view name) {
    // name.strip().lstrip("0123456789")
    std::string n = core::pythonStrip(name);
    n.erase(0, std::min(n.find_first_not_of("0123456789"), n.size()));
    std::string upper = n;
    for (char& c : upper) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    for (const char* two : {"MG", "NA", "CL", "ZN", "CA"}) {
        if (upper == two) {  // startswith(two) and the same length
            return {two[0], static_cast<char>(two[1] - 'A' + 'a')};
        }
    }
    if (n.empty()) {
        return "X";
    }
    // name[:1].upper(): the first character (code point)
    const auto lead = static_cast<unsigned char>(n[0]);
    std::size_t len = 1;
    if (lead >= 0xF0) {
        len = 4;
    } else if (lead >= 0xE0) {
        len = 3;
    } else if (lead >= 0xC0) {
        len = 2;
    }
    return upper.substr(0, std::min(len, upper.size()));
}

namespace {

using core::NumArray;

core::Array stringArray(const std::vector<std::string>& items) {
    return core::Array(items.begin(), items.end());
}

constexpr double kLinkMax = 9.0;  // LINK_MAX: max trace-point distance (Å) of a drawn backbone link

core::AtomGroup residueAtomGroup(const core::AtomGroup& like, std::int32_t residue) {
    const auto atoms = like.topology().residueAtoms(residue);
    return core::AtomGroup(like.sharedTopology(), std::vector<std::int32_t>(atoms.begin(), atoms.end()));
}

Vec3 position(const core::Frame& frame, std::int32_t atom) {
    const auto k = static_cast<std::size_t>(atom) * 3;
    return {frame.positions[k], frame.positions[k + 1], frame.positions[k + 2]};
}

// The positions of the atoms, widened (atoms.positions.astype(np.float64)).
std::vector<Vec3> positions(const core::Frame& frame, const core::AtomGroup& atoms) {
    std::vector<Vec3> out;
    out.reserve(atoms.size());
    for (const std::int32_t a : atoms.indices()) {
        out.push_back(position(frame, a));
    }
    return out;
}

// The float32 positions of the atoms, x, y, z interleaved (atoms.positions).
std::vector<float> positions32(const core::Frame& frame, const core::AtomGroup& atoms) {
    std::vector<float> out;
    out.reserve(atoms.size() * 3);
    for (const std::int32_t a : atoms.indices()) {
        const auto k = static_cast<std::size_t>(a) * 3;
        out.insert(out.end(), frame.positions.begin() + static_cast<std::ptrdiff_t>(k),
                   frame.positions.begin() + static_cast<std::ptrdiff_t>(k + 3));
    }
    return out;
}

const core::Frame& currentFrame(const AnalysisContext& context) {
    const core::Frame* frame = context.system().currentFrame();
    if (!frame) {
        throw core::ValueError("This Universe has no coordinates");
    }
    return *frame;
}

// (x - mc) @ r.T + rc, row by row
void superpose(std::vector<Vec3>& xyz, const KabschResult& fit) {
    for (Vec3& x : xyz) {
        const Vec3 d{x[0] - fit.mobileCenter[0], x[1] - fit.mobileCenter[1], x[2] - fit.mobileCenter[2]};
        Vec3 y{};
        for (std::size_t j = 0; j < 3; ++j) {
            const auto& r = fit.rotation[j];
            y[j] = ((d[0] * r[0] + d[1] * r[1]) + d[2] * r[2]) + fit.referenceCenter[j];
        }
        x = y;
    }
}

void store(std::vector<float>& out, std::size_t k, const std::vector<Vec3>& xyz) {
    std::size_t at = k * xyz.size() * 3;
    for (const Vec3& x : xyz) {
        for (const double v : x) {
            out[at++] = static_cast<float>(v);
        }
    }
}

// residue_group_membership: [n_groups, n_residues], row-major; groups whose selection fails are
// left out (all false).
std::vector<bool> groupMembership(const AnalysisContext& context, const std::vector<std::int32_t>& residues) {
    const std::vector<core::ResidueGroup>& groups = context.project().groups;
    std::vector<bool> out(groups.size() * residues.size(), false);
    for (std::size_t g = 0; g < groups.size(); ++g) {
        core::AtomGroup atoms;
        try {
            atoms = context.select("group:" + groups[g].name, true);
        } catch (const core::SelectionError&) {
            continue;
        }
        const std::vector<std::int32_t> in = atoms.residues();
        for (std::size_t r = 0; r < residues.size(); ++r) {
            out[g * residues.size() + r] = std::binary_search(in.begin(), in.end(), residues[r]);
        }
    }
    return out;
}

// residue_group_index: the first group containing each residue, -1 for none.
std::vector<std::int64_t> groupIndex(const std::vector<bool>& member, std::size_t nGroups, std::size_t nResidues) {
    std::vector<std::int64_t> out(nResidues, -1);
    for (std::size_t g = nGroups; g-- > 0;) {
        for (std::size_t r = 0; r < nResidues; ++r) {
            if (member[g * nResidues + r]) {
                out[r] = static_cast<std::int64_t>(g);
            }
        }
    }
    return out;
}

// TracePoints: one point per residue, the first atom of trace_atoms present (C4', CA by default),
// else the centroid of the residue's heavy atoms (of all its atoms without heavy atoms).
class TracePoints {
public:
    TracePoints(const core::AtomGroup& like, const std::vector<std::int32_t>& residues, std::string_view traceAtoms) {
        std::vector<std::string> priority = core::splitCsv(traceAtoms);
        if (priority.empty()) {
            priority = {"C4'", "CA"};
        }
        const core::Topology& t = like.topology();
        for (std::size_t row = 0; row < residues.size(); ++row) {
            // {a.name: a.ix}: the last atom of a name wins
            std::unordered_map<std::string, std::int32_t> names;
            for (const std::int32_t a : t.residueAtoms(residues[row])) {
                names[t.name(a)] = a;
            }
            std::int64_t chosen = -1;
            for (const std::string& n : priority) {
                if (const auto it = names.find(n); it != names.end()) {
                    chosen = it->second;
                    break;
                }
            }
            if (chosen < 0) {
                const core::AtomGroup all = residueAtomGroup(like, residues[row]);
                const core::AtomGroup hv = core::heavy(all);
                m_centroidRows.push_back(row);
                m_centroidAtoms.push_back(hv.empty() ? all.indices() : hv.indices());
            }
            m_pointAtoms.push_back(chosen);
            m_isPoint.push_back(chosen >= 0);
        }
    }

    std::size_t size() const { return m_pointAtoms.size(); }
    const std::vector<bool>& isPoint() const { return m_isPoint; }

    std::vector<Vec3> positions(const core::Frame& frame) const {
        std::vector<Vec3> xyz(size(), Vec3{0.0, 0.0, 0.0});
        for (std::size_t row = 0; row < size(); ++row) {
            if (m_pointAtoms[row] >= 0) {
                xyz[row] = position(frame, static_cast<std::int32_t>(m_pointAtoms[row]));
            }
        }
        for (std::size_t c = 0; c < m_centroidRows.size(); ++c) {
            // ag.positions.mean(axis=0) on float32: rows accumulated in turn in float32, then divided
            // by the count (numpy divides by an intp: in float64, then rounded to float32, which
            // equals the float32 quotient)
            const std::vector<std::int32_t>& atoms = m_centroidAtoms[c];
            std::array<float, 3> sum{0.0f, 0.0f, 0.0f};
            for (const std::int32_t a : atoms) {
                const auto k = static_cast<std::size_t>(a) * 3;
                for (std::size_t d = 0; d < 3; ++d) {
                    sum[d] = sum[d] + frame.positions[k + d];
                }
            }
            Vec3& out = xyz[m_centroidRows[c]];
            for (std::size_t d = 0; d < 3; ++d) {
                out[d] = static_cast<float>(static_cast<double>(sum[d]) / static_cast<double>(atoms.size()));
            }
        }
        return xyz;
    }

private:
    std::vector<std::int64_t> m_pointAtoms;
    std::vector<bool> m_isPoint;
    std::vector<std::size_t> m_centroidRows;
    std::vector<std::vector<std::int32_t>> m_centroidAtoms;
};

// atom_bonds(atoms, 1.95): bonds among `atoms` as local index pairs [m, 2]; from the topology when
// it has bonds among them (sorted, unique), else guessed from the current positions: pairs within
// 1.95 Å in the same or adjacent residues (by residue index), never H-H, and with an H only within
// 1.25 Å. The guessed rows are sorted by (i, j) (MDATK: kd-tree order; see the file comment).
std::vector<std::int64_t> atomBonds(const core::AtomGroup& atoms, const core::Frame& frame) {
    constexpr double maxLen = 1.95;
    const core::Topology& t = atoms.topology();
    std::unordered_map<std::int32_t, std::int64_t> index;
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        index[atoms[k]] = static_cast<std::int64_t>(k);
    }
    if (t.bonds) {
        std::set<std::pair<std::int64_t, std::int64_t>> pairs;
        for (const core::Bond& b : atoms.bonds()) {
            const auto i = index.find(b[0]);
            const auto j = index.find(b[1]);
            if (i != index.end() && j != index.end()) {
                pairs.insert(std::minmax(i->second, j->second));
            }
        }
        if (!pairs.empty()) {
            std::vector<std::int64_t> out;
            for (const auto& [i, j] : pairs) {
                out.push_back(i);
                out.push_back(j);
            }
            return out;
        }
    }
    if (atoms.size() < 2) {
        return {};
    }
    std::vector<double> xyz;
    xyz.reserve(atoms.size() * 3);
    for (const Vec3& x : positions(frame, atoms)) {
        xyz.insert(xyz.end(), x.begin(), x.end());
    }
    const PairList found = pairsWithin(xyz, xyz, maxLen, std::nullopt);
    std::vector<bool> h;
    h.reserve(atoms.size());
    for (const std::int32_t a : atoms.indices()) {
        h.push_back(structureGuessElement(t.name(a)) == "H");
    }
    std::vector<std::int64_t> out;
    for (std::size_t p = 0; p < found.size(); ++p) {
        const std::int64_t i = found.i[p];
        const std::int64_t j = found.j[p];
        if (i >= j) {
            continue;
        }
        const auto ui = static_cast<std::size_t>(i);
        const auto uj = static_cast<std::size_t>(j);
        const std::int64_t ri = t.atomResidue[static_cast<std::size_t>(atoms[ui])];
        const std::int64_t rj = t.atomResidue[static_cast<std::size_t>(atoms[uj])];
        const bool near = std::abs(ri - rj) <= 1;
        const bool hh = h[ui] && h[uj];
        if (near && !hh && (found.d[p] <= 1.25 || !(h[ui] || h[uj]))) {
            out.push_back(i);
            out.push_back(j);
        }
    }
    return out;
}

std::vector<std::int64_t> residueRowsOf(const std::vector<std::int32_t>& residues, const core::AtomGroup& atoms) {
    // np.searchsorted(residues.ix, atoms.resindices)
    std::vector<std::int64_t> out;
    out.reserve(atoms.size());
    for (const std::int32_t a : atoms.indices()) {
        const std::int32_t r = atoms.topology().atomResidue[static_cast<std::size_t>(a)];
        out.push_back(static_cast<std::int64_t>(std::lower_bound(residues.begin(), residues.end(), r) -
                                                residues.begin()));
    }
    return out;
}

std::vector<std::int64_t> widen(const std::vector<std::int32_t>& values) {
    return {values.begin(), values.end()};
}

std::vector<std::int64_t> resids(const core::Topology& t, const std::vector<std::int32_t>& residues) {
    std::vector<std::int64_t> out;
    out.reserve(residues.size());
    for (const std::int32_t r : residues) {
        out.push_back(t.resids[static_cast<std::size_t>(r)]);
    }
    return out;
}

// ---------------------------------------------------------------- trace

class Trace final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::Project& p = project();
        StructureResidues structure = structureResidues(context());
        if (structure.residues.empty()) {
            throw core::SelectionError("structural view: no residues selected");
        }
        m_residues = std::move(structure.residues);
        const core::AtomGroup all = context().system().atoms();
        m_tracer.emplace(all, m_residues, p.visual.trace_atoms);
        m_ions = p.visual.show_ions ? context().select(p.analyses.mg.mg_selection, true)
                                    : core::AtomGroup(all.sharedTopology(), {});
        m_align = p.visual.align_frames;
        const std::size_t n = frameCount();
        m_coords.assign(n * m_residues.size() * 3, 0.0f);
        m_ionCoords.assign(n * m_ions.size() * 3, 0.0f);
        m_reference.reset();
        const core::Topology& t = context().topology();
        m_chain = core::chainIds(t, m_residues);
        m_groupMember = groupMembership(context(), m_residues);
        m_groupIndex = groupIndex(m_groupMember, p.groups.size(), m_residues.size());
        m_links.clear();
    }

    // The first analysed frame is the alignment reference and decides the links.
    void prepare(const core::Frame& first) override {
        m_reference = m_tracer->positions(first);
        links(*m_reference);
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        std::vector<Vec3> xyz = m_tracer->positions(frame);
        std::vector<Vec3> ions = positions(frame, m_ions);
        if (k != 0 && m_align && m_residues.size() >= 3) {  // the reference frame stays as it is
            const KabschResult fit = kabsch(xyz, *m_reference);
            superpose(xyz, fit);
            superpose(ions, fit);
        }
        store(m_coords, k, xyz);
        store(m_ionCoords, k, ions);
    }

    core::AnalysisResult finalize() override {
        const core::Topology& t = context().topology();
        std::vector<std::string> ionLabels;
        for (const std::int32_t a : m_ions.indices()) {
            const std::int32_t r = t.atomResidue[static_cast<std::size_t>(a)];
            ionLabels.push_back(t.resname(r) + std::to_string(t.resids[static_cast<std::size_t>(r)]));
        }
        const std::vector<std::string> groups = core::groupNames(project());
        core::AnalysisResult res = newResult({}, {{"residue_labels", stringArray(core::residueLabels(t, m_residues))},
                                                  {"ion_labels", stringArray(ionLabels)},
                                                  {"group_names", stringArray(groups)}});
        const std::size_t n = frameCount();
        const std::size_t nRes = m_residues.size();
        std::vector<std::string> resnames;
        for (const std::int32_t r : m_residues) {
            resnames.push_back(t.resname(r));
        }
        res.arrays.emplace_back("coords", NumArray(m_coords, {n, nRes, 3}));
        res.arrays.emplace_back("ion_coords", NumArray(m_ionCoords, {n, m_ions.size(), 3}));
        res.arrays.emplace_back("chain", NumArray(m_chain));
        res.arrays.emplace_back("group_index", NumArray(m_groupIndex));
        res.arrays.emplace_back("group_member", NumArray(m_groupMember, {project().groups.size(), nRes}));
        res.arrays.emplace_back("is_point", NumArray(m_tracer->isPoint()));
        res.arrays.emplace_back("links", NumArray(m_links, {m_links.size() / 2, 2}));
        res.arrays.emplace_back("resindex", NumArray(widen(m_residues)));
        res.arrays.emplace_back("resid", NumArray(resids(t, m_residues)));
        res.arrays.emplace_back("resname", NumArray::strings(resnames, 8));
        res.arrays.emplace_back("chain_name", NumArray(core::chainNames(t, m_residues)));
        return res;
    }

private:
    // Links between consecutive residues (by residue index) of one chain whose trace points are
    // atoms no farther apart than LINK_MAX, on the first frame.
    void links(const std::vector<Vec3>& xyz) {
        const std::vector<bool>& isPoint = m_tracer->isPoint();
        for (std::size_t i = 0; i + 1 < xyz.size(); ++i) {
            const std::size_t j = i + 1;
            const double dx = xyz[j][0] - xyz[i][0];
            const double dy = xyz[j][1] - xyz[i][1];
            const double dz = xyz[j][2] - xyz[i][2];
            const bool close = std::sqrt((dx * dx + dy * dy) + dz * dz) <= kLinkMax;
            if (m_chain[i] == m_chain[j] && m_residues[j] - m_residues[i] == 1 && isPoint[i] && isPoint[j] && close) {
                m_links.push_back(static_cast<std::int64_t>(i));
                m_links.push_back(static_cast<std::int64_t>(j));
            }
        }
    }

    std::vector<std::int32_t> m_residues;
    std::optional<TracePoints> m_tracer;
    core::AtomGroup m_ions;
    bool m_align = true;
    std::vector<float> m_coords;
    std::vector<float> m_ionCoords;
    std::optional<std::vector<Vec3>> m_reference;
    std::vector<std::int64_t> m_chain;
    std::vector<bool> m_groupMember;
    std::vector<std::int64_t> m_groupIndex;
    std::vector<std::int64_t> m_links;
};

// ---------------------------------------------------------------- proximity

class Proximity final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::Project& p = project();
        StructureResidues structure = structureResidues(context());
        m_residues = std::move(structure.residues);
        const std::string text = core::resolveSelection(p.visual.proximity_set, p.groups, &context().topology());
        m_set = core::heavy(context().select("sel:" + text));
        const std::vector<std::int32_t> setResidues = m_set.residues();
        m_member.clear();
        for (const std::int32_t r : m_residues) {
            m_member.push_back(std::binary_search(setResidues.begin(), setResidues.end(), r));
        }
        const core::AtomGroup heavyAtoms = core::heavy(structure.atoms);
        std::vector<bool> keep;
        keep.reserve(heavyAtoms.size());
        for (const std::int32_t a : heavyAtoms.indices()) {
            const std::int32_t r = context().topology().atomResidue[static_cast<std::size_t>(a)];
            const auto it = std::lower_bound(m_residues.begin(), m_residues.end(), r);
            keep.push_back(!m_member[static_cast<std::size_t>(it - m_residues.begin())]);
        }
        m_struct = heavyAtoms.masked(keep);
        m_row = residueRowsOf(m_residues, m_struct);
        m_radius = p.analyses.proximity_search_radius;
        m_dist.assign(m_residues.size() * frameCount(), std::numeric_limits<float>::infinity());
    }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        if (m_struct.empty()) {
            return;
        }
        const std::size_t n = frameCount();
        const PairList found = pairsWithin(std::span<const float>(positions32(frame, m_struct)),
                                           std::span<const float>(positions32(frame, m_set)), m_radius,
                                           usableBox(context().box(frame)));
        if (!found.empty()) {
            std::vector<double> best(m_residues.size(), std::numeric_limits<double>::infinity());
            for (std::size_t p = 0; p < found.size(); ++p) {
                double& b = best[static_cast<std::size_t>(m_row[static_cast<std::size_t>(found.i[p])])];
                b = std::min(b, found.d[p]);  // np.minimum.at
            }
            for (std::size_t r = 0; r < best.size(); ++r) {
                m_dist[r * n + k] = static_cast<float>(best[r]);
            }
        }
        for (std::size_t r = 0; r < m_member.size(); ++r) {
            if (m_member[r]) {
                m_dist[r * n + k] = 0.0f;
            }
        }
    }

    core::AnalysisResult finalize() override {
        const core::Topology& t = context().topology();
        core::AnalysisResult res =
            newResult({}, {{"residue_labels", stringArray(core::residueLabels(t, m_residues))}, {"radius", m_radius}});
        res.arrays.emplace_back("min_distance", NumArray(m_dist, {m_residues.size(), frameCount()}));
        res.arrays.emplace_back("member", NumArray(m_member));
        res.arrays.emplace_back("resindex", NumArray(widen(m_residues)));
        res.addMatrix("min_distance", "min_distance", "residue_labels", "Distance to proximity set", "value",
                      {{"unit", "Å"}});
        return res;
    }

private:
    std::vector<std::int32_t> m_residues;
    core::AtomGroup m_set;
    std::vector<bool> m_member;
    core::AtomGroup m_struct;
    std::vector<std::int64_t> m_row;
    double m_radius = 0.0;
    std::vector<float> m_dist;  // [residues, frames]
};

// ---------------------------------------------------------------- atoms

class Atoms final : public Analysis {
public:
    using Analysis::Analysis;

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        const core::Project& p = project();
        StructureResidues structure = structureResidues(context());
        if (structure.residues.empty()) {
            throw core::SelectionError("structural view: no residues selected");
        }
        core::AtomGroup atoms = p.visual.rep_hydrogens ? structure.atoms : core::heavy(structure.atoms);
        std::vector<std::int32_t> order = atoms.indices();
        std::stable_sort(order.begin(), order.end());  // np.argsort(atoms.ix, kind="stable")
        m_atoms = core::AtomGroup(atoms.sharedTopology(), std::move(order));
        m_residues = std::move(structure.residues);
        m_tracer.emplace(context().system().atoms(), m_residues, p.visual.trace_atoms);
        m_align = p.visual.align_frames;
        m_reference.reset();
        m_coords.assign(frameCount() * m_atoms.size() * 3, 0.0f);
        m_row = residueRowsOf(m_residues, m_atoms);
        m_bonds = atomBonds(m_atoms, currentFrame(context()));  // the frame current at setup
    }

    // The first analysed frame is the alignment reference.
    void prepare(const core::Frame& first) override { m_reference = m_tracer->positions(first); }

    bool concurrentProcess() const override { return true; }

    void process(std::size_t k, const core::Frame& frame) override {
        std::vector<Vec3> xyz = positions(frame, m_atoms);
        if (k != 0 && m_align && m_tracer->size() >= 3) {  // the reference frame stays as it is
            superpose(xyz, kabsch(m_tracer->positions(frame), *m_reference));
        }
        store(m_coords, k, xyz);
    }

    core::AnalysisResult finalize() override {
        const core::Topology& t = context().topology();
        std::vector<std::string> names;
        std::vector<std::string> resnames;
        std::vector<std::int64_t> atomResids;
        for (const std::int32_t a : m_atoms.indices()) {
            const std::int32_t r = t.atomResidue[static_cast<std::size_t>(a)];
            names.push_back(t.name(a));
            resnames.push_back(t.resname(r));
            atomResids.push_back(t.resids[static_cast<std::size_t>(r)]);
        }
        core::AnalysisResult res = newResult({}, {{"atom_names", stringArray(names)},
                                                  {"elements", stringArray(core::atomElements(m_atoms))},
                                                  {"resnames", stringArray(resnames)}});
        res.arrays.emplace_back("coords", NumArray(m_coords, {frameCount(), m_atoms.size(), 3}));
        res.arrays.emplace_back("row", NumArray(m_row));
        res.arrays.emplace_back("bonds", NumArray(m_bonds, {m_bonds.size() / 2, 2}));
        res.arrays.emplace_back("resid", NumArray(std::move(atomResids)));
        std::vector<std::int64_t> atomIndex(m_atoms.indices().begin(), m_atoms.indices().end());
        res.arrays.emplace_back("atom_index", NumArray(std::move(atomIndex)));
        res.arrays.emplace_back("resindex", NumArray(widen(m_residues)));
        return res;
    }

private:
    core::AtomGroup m_atoms;
    std::vector<std::int32_t> m_residues;
    std::optional<TracePoints> m_tracer;
    bool m_align = true;
    std::optional<std::vector<Vec3>> m_reference;
    std::vector<float> m_coords;
    std::vector<std::int64_t> m_row;
    std::vector<std::int64_t> m_bonds;  // [m, 2] row-major
};

// ---------------------------------------------------------------- descriptors

class TraceDescriptor final : public AnalysisDescriptor {
public:
    // version 4: v2 residue numbers for tags; v3 group membership; v4 chain and residue names
    TraceDescriptor() : AnalysisDescriptor({"trace", "Backbone trace", 4, true, 0, {}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::VisualParams& v = project.visual;
        return {{"groups", stringArray(structureGroups(project))},
                {"trace_atoms", v.trace_atoms},
                {"align", v.align_frames},
                {"show_ions", v.show_ions},
                {"ions", v.show_ions ? project.analyses.mg.mg_selection : std::string()},
                {"group_order", stringArray(core::groupNames(project))}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<Trace>(*this, context, std::move(instance));
    }
};

class ProximityDescriptor final : public AnalysisDescriptor {
public:
    ProximityDescriptor() : AnalysisDescriptor({"proximity", "Proximity to residue set", 1, true, 0, {}}) {}

    ExpectedOutputs expectedOutputs(const core::Project&) const override { return {{}, {"min_distance"}, {}}; }

    core::Table paramsFor(const core::Project& project) const override {
        return {{"groups", stringArray(structureGroups(project))},
                {"set", project.visual.proximity_set},
                {"radius", project.analyses.proximity_search_radius},
                {"ions", project.analyses.mg.mg_selection}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.visual.proximity_set);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<Proximity>(*this, context, std::move(instance));
    }
};

class AtomsDescriptor final : public AnalysisDescriptor {
public:
    // version 2: elements guessed per residue (protein CA is carbon, not calcium)
    AtomsDescriptor() : AnalysisDescriptor({"atoms", "Atomic coordinates (structure view)", 2, true, 0, {}}) {}

    core::Table paramsFor(const core::Project& project) const override {
        const core::VisualParams& v = project.visual;
        return {{"groups", stringArray(structureGroups(project))},
                {"trace_atoms", v.trace_atoms},
                {"align", v.align_frames},
                {"hydrogens", v.rep_hydrogens},
                {"ions", project.analyses.mg.mg_selection},
                {"group_order", stringArray(core::groupNames(project))}};
    }

    std::vector<std::string> selectionTexts(const core::Project& project) const override {
        std::vector<std::string> out = structureGroups(project);
        out.push_back(project.analyses.mg.mg_selection);
        return out;
    }

    std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance) const override {
        return std::make_unique<Atoms>(*this, context, std::move(instance));
    }
};

}  // namespace

const AnalysisDescriptor& traceAnalysis() {
    static const TraceDescriptor descriptor;
    return descriptor;
}

const AnalysisDescriptor& proximityAnalysis() {
    static const ProximityDescriptor descriptor;
    return descriptor;
}

const AnalysisDescriptor& atomsAnalysis() {
    static const AtomsDescriptor descriptor;
    return descriptor;
}

}  // namespace mdbridge::analysis
