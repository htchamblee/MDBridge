#include "mdbridge/core/topology.hpp"

#include <algorithm>
#include <numeric>
#include <set>
#include <stdexcept>

namespace mdbridge::core {

namespace {

template <class T>
void checkSize(const std::optional<std::vector<T>>& values, std::size_t expected, const char* what) {
    if (values && values->size() != expected) {
        throw std::invalid_argument(std::string("topology: ") + what + " has " + std::to_string(values->size()) +
                                    " entries, expected " + std::to_string(expected));
    }
}

template <class T>
void checkSize(const std::vector<T>& values, std::size_t expected, const char* what) {
    if (values.size() != expected) {
        throw std::invalid_argument(std::string("topology: ") + what + " has " + std::to_string(values.size()) +
                                    " entries, expected " + std::to_string(expected));
    }
}

void checkStrings(const std::vector<std::uint32_t>& ids, const StringPool& strings, const char* what) {
    for (const std::uint32_t id : ids) {
        if (id >= strings.size()) {
            throw std::invalid_argument(std::string("topology: ") + what + " refers to an unknown string");
        }
    }
}

void checkStrings(const std::optional<std::vector<std::uint32_t>>& ids, const StringPool& strings, const char* what) {
    if (ids) {
        checkStrings(*ids, strings, what);
    }
}

// first atom of each residue (np.unique(residues.atoms.resindices, return_index=True)): the
// residues in ascending order, each once
std::vector<std::int32_t> firstAtoms(const Topology& topology, std::span<const std::int32_t> residues) {
    std::set<std::int32_t> unique(residues.begin(), residues.end());
    std::vector<std::int32_t> atoms;
    atoms.reserve(unique.size());
    for (const std::int32_t r : unique) {
        const auto members = topology.residueAtoms(r);
        if (!members.empty()) {
            atoms.push_back(members.front());
        }
    }
    return atoms;
}

// np.unique(values, return_inverse=True)[1]: the rank of each value among the sorted distinct values
template <class T>
std::vector<std::int64_t> ranks(const std::vector<T>& values) {
    std::vector<T> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::vector<std::int64_t> out;
    out.reserve(values.size());
    for (const T& v : values) {
        out.push_back(std::lower_bound(sorted.begin(), sorted.end(), v) - sorted.begin());
    }
    return out;
}

template <class T>
std::size_t distinct(const std::vector<T>& values) {
    return std::set<T>(values.begin(), values.end()).size();
}

std::string pythonStrip(const std::string& text) {
    // chain IDs and segids are ASCII in practice; str(c).strip() of the C++ strings
    const auto first = text.find_first_not_of(" \t\n\r\v\f");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = text.find_last_not_of(" \t\n\r\v\f");
    return text.substr(first, last - first + 1);
}

}  // namespace

StringPool::StringPool() { intern(""); }

std::uint32_t StringPool::intern(std::string_view text) {
    if (const auto it = m_ids.find(text); it != m_ids.end()) {
        return it->second;
    }
    const auto id = static_cast<std::uint32_t>(m_strings.size());
    m_strings.emplace_back(text);
    m_ids.emplace(std::string(text), id);
    return id;
}

std::optional<std::uint32_t> StringPool::find(std::string_view text) const {
    if (const auto it = m_ids.find(text); it != m_ids.end()) {
        return it->second;
    }
    return std::nullopt;
}

void Topology::finish() {
    const std::size_t n = names.size();
    if (types.empty() && n > 0) {
        types.assign(n, 0);
    }
    checkSize(types, n, "types");
    checkSize(atomResidue, n, "atomResidue");
    checkSize(elements, n, "elements");
    checkSize(chainIDs, n, "chainIDs");
    checkSize(masses, n, "masses");
    checkSize(charges, n, "charges");
    checkSize(ids, n, "ids");
    checkSize(altLocs, n, "altLocs");
    checkSize(recordTypes, n, "recordTypes");
    checkSize(occupancies, n, "occupancies");
    checkSize(tempFactors, n, "tempFactors");
    checkSize(formalCharges, n, "formalCharges");

    const std::size_t nr = resids.size();
    if (resnums.empty() && nr > 0) {
        resnums = resids;
    }
    checkSize(resnums, nr, "resnums");
    checkSize(resnames, nr, "resnames");
    checkSize(icodes, nr, "icodes");
    checkSize(residueSegment, nr, "residueSegment");

    checkStrings(names, strings, "names");
    checkStrings(types, strings, "types");
    checkStrings(elements, strings, "elements");
    checkStrings(chainIDs, strings, "chainIDs");
    checkStrings(altLocs, strings, "altLocs");
    checkStrings(recordTypes, strings, "recordTypes");
    checkStrings(resnames, strings, "resnames");
    checkStrings(icodes, strings, "icodes");
    checkStrings(segids, strings, "segids");

    for (const std::int32_t r : atomResidue) {
        if (r < 0 || static_cast<std::size_t>(r) >= nr) {
            throw std::invalid_argument("topology: an atom refers to a residue that does not exist");
        }
    }
    for (const std::int32_t s : residueSegment) {
        if (s < 0 || static_cast<std::size_t>(s) >= segids.size()) {
            throw std::invalid_argument("topology: a residue refers to a segment that does not exist");
        }
    }
    if (bonds) {
        for (Bond& bond : *bonds) {
            if (bond[0] > bond[1]) {
                std::swap(bond[0], bond[1]);
            }
            if (bond[0] < 0 || static_cast<std::size_t>(bond[1]) >= n || bond[0] == bond[1]) {
                throw std::invalid_argument("topology: a bond refers to atoms that do not exist");
            }
        }
        std::sort(bonds->begin(), bonds->end());
        bonds->erase(std::unique(bonds->begin(), bonds->end()), bonds->end());
    }

    m_residueAtomOffsets.assign(nr + 1, 0);
    for (const std::int32_t r : atomResidue) {
        ++m_residueAtomOffsets[static_cast<std::size_t>(r) + 1];
    }
    std::partial_sum(m_residueAtomOffsets.begin(), m_residueAtomOffsets.end(), m_residueAtomOffsets.begin());
    m_residueAtoms.assign(n, 0);
    std::vector<std::int32_t> fill(m_residueAtomOffsets.begin(), m_residueAtomOffsets.end() - 1);
    for (std::size_t atom = 0; atom < n; ++atom) {
        m_residueAtoms[static_cast<std::size_t>(fill[static_cast<std::size_t>(atomResidue[atom])]++)] =
            static_cast<std::int32_t>(atom);
    }
}

std::span<const std::int32_t> Topology::residueAtoms(std::int32_t residue) const {
    const auto r = static_cast<std::size_t>(residue);
    if (r + 1 >= m_residueAtomOffsets.size()) {
        return {};
    }
    return std::span<const std::int32_t>(m_residueAtoms)
        .subspan(static_cast<std::size_t>(m_residueAtomOffsets[r]),
                 static_cast<std::size_t>(m_residueAtomOffsets[r + 1] - m_residueAtomOffsets[r]));
}

std::vector<std::int32_t> fragmentIndices(const Topology& topology) {
    if (!topology.bonds) {
        throw NoDataError("This Universe does not contain bonds information");
    }
    const std::size_t n = topology.atomCount();
    std::vector<std::int32_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&parent](std::int32_t a) {
        while (parent[static_cast<std::size_t>(a)] != a) {
            parent[static_cast<std::size_t>(a)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
            a = parent[static_cast<std::size_t>(a)];
        }
        return a;
    };
    for (const Bond& bond : *topology.bonds) {
        const std::int32_t a = root(bond[0]);
        const std::int32_t b = root(bond[1]);
        if (a != b) {
            parent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);  // the root is the lowest atom
        }
    }
    std::vector<std::int32_t> label(n, -1);
    std::vector<std::int32_t> out(n);
    std::int32_t next = 0;
    for (std::size_t atom = 0; atom < n; ++atom) {
        const auto r = static_cast<std::size_t>(root(static_cast<std::int32_t>(atom)));
        if (label[r] < 0) {
            label[r] = next++;
        }
        out[atom] = label[r];
    }
    return out;
}

AtomGroup::AtomGroup(std::shared_ptr<const Topology> topology, std::vector<std::int32_t> indices)
    : m_topology(std::move(topology)), m_indices(std::move(indices)) {}

AtomGroup::AtomGroup(std::shared_ptr<const Topology> topology) : m_topology(std::move(topology)) {
    if (m_topology) {
        m_indices.resize(m_topology->atomCount());
        std::iota(m_indices.begin(), m_indices.end(), 0);
    }
}

AtomGroup AtomGroup::masked(const std::vector<bool>& mask) const {
    if (mask.size() != m_indices.size()) {
        throw std::invalid_argument("atom group: the mask has " + std::to_string(mask.size()) + " entries for " +
                                    std::to_string(m_indices.size()) + " atoms");
    }
    std::vector<std::int32_t> out;
    for (std::size_t i = 0; i < m_indices.size(); ++i) {
        if (mask[i]) {
            out.push_back(m_indices[i]);
        }
    }
    return AtomGroup(m_topology, std::move(out));
}

AtomGroup AtomGroup::at(std::span<const std::int64_t> positions) const {
    std::vector<std::int32_t> out;
    out.reserve(positions.size());
    const auto n = static_cast<std::int64_t>(m_indices.size());
    for (std::int64_t p : positions) {
        if (p < 0) {
            p += n;
        }
        if (p < 0 || p >= n) {
            throw std::out_of_range("index " + std::to_string(p) + " is out of bounds for size " + std::to_string(n));
        }
        out.push_back(m_indices[static_cast<std::size_t>(p)]);
    }
    return AtomGroup(m_topology, std::move(out));
}

AtomGroup AtomGroup::unique() const {
    std::vector<std::int32_t> out = m_indices;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return AtomGroup(m_topology, std::move(out));
}

std::vector<std::int32_t> AtomGroup::residues() const {
    std::vector<std::int32_t> out;
    out.reserve(m_indices.size());
    for (const std::int32_t atom : m_indices) {
        out.push_back(m_topology->atomResidue[static_cast<std::size_t>(atom)]);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<Bond> AtomGroup::bonds() const {
    if (!m_topology->bonds) {
        throw NoDataError("This Universe does not contain bonds information");
    }
    std::vector<bool> inside(m_topology->atomCount(), false);
    for (const std::int32_t atom : m_indices) {
        inside[static_cast<std::size_t>(atom)] = true;
    }
    std::vector<Bond> out;
    for (const Bond& bond : *m_topology->bonds) {
        if (inside[static_cast<std::size_t>(bond[0])] || inside[static_cast<std::size_t>(bond[1])]) {
            out.push_back(bond);
        }
    }
    return out;
}

AtomGroup operator|(const AtomGroup& a, const AtomGroup& b) {
    if (a.sharedTopology() != b.sharedTopology()) {
        throw std::invalid_argument("atom groups of different topologies cannot be combined");
    }
    std::vector<std::int32_t> out = a.indices();
    out.insert(out.end(), b.indices().begin(), b.indices().end());
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return AtomGroup(a.sharedTopology(), std::move(out));
}

AtomGroup residueAtoms(const std::shared_ptr<const Topology>& topology, std::span<const std::int32_t> residues) {
    std::vector<std::int32_t> out;
    for (const std::int32_t r : residues) {
        const auto atoms = topology->residueAtoms(r);
        out.insert(out.end(), atoms.begin(), atoms.end());
    }
    return AtomGroup(topology, std::move(out));
}

std::vector<std::int64_t> chainIds(const Topology& topology, std::span<const std::int32_t> residues) {
    const std::size_t n = residues.size();
    if (n == 0) {
        return {};
    }
    const std::vector<std::int32_t> first = firstAtoms(topology, residues);
    if (topology.chainIDs) {
        std::vector<std::string> cid;
        for (const std::int32_t atom : first) {
            cid.push_back(topology.strings[(*topology.chainIDs)[static_cast<std::size_t>(atom)]]);
        }
        if (distinct(cid) > 1) {
            return ranks(cid);
        }
    }
    std::vector<std::string> segs;
    for (const std::int32_t r : residues) {
        segs.push_back(topology.segid(topology.residueSegment[static_cast<std::size_t>(r)]));
    }
    if (distinct(segs) > 1) {
        return ranks(segs);
    }
    if (!topology.bonds) {
        return std::vector<std::int64_t>(n, 0);
    }
    const std::vector<std::int32_t> fragments = fragmentIndices(topology);
    std::vector<std::int32_t> frag;
    for (const std::int32_t atom : first) {
        frag.push_back(fragments[static_cast<std::size_t>(atom)]);
    }
    return ranks(frag);
}

std::vector<std::string> chainNames(const Topology& topology, std::span<const std::int32_t> residues) {
    const std::size_t n = residues.size();
    if (n == 0) {
        return {};
    }
    std::vector<std::string> cid;
    if (topology.chainIDs) {
        for (const std::int32_t atom : firstAtoms(topology, residues)) {
            cid.push_back(pythonStrip(topology.strings[(*topology.chainIDs)[static_cast<std::size_t>(atom)]]));
        }
    } else {
        cid.assign(n, "");
    }
    std::vector<std::string> segs;
    for (const std::int32_t r : residues) {
        segs.push_back(pythonStrip(topology.segid(topology.residueSegment[static_cast<std::size_t>(r)])));
    }
    if (distinct(cid) > 1) {
        return cid;
    }
    if (distinct(segs) > 1) {
        return segs;
    }
    return cid;
}

MemoryTrajectory::MemoryTrajectory(std::vector<Frame> frames, double dt, bool keepTimes)
    : m_frames(std::move(frames)), m_dt(dt) {
    for (std::size_t i = 0; i < m_frames.size(); ++i) {
        m_frames[i].index = static_cast<std::int64_t>(i);
        if (!keepTimes) {
            m_frames[i].time = static_cast<double>(i) * dt;
        }
    }
}

std::int64_t MemoryTrajectory::frameCount() const { return static_cast<std::int64_t>(m_frames.size()); }

const Frame& MemoryTrajectory::seek(std::int64_t index) {
    const auto n = frameCount();
    if (index < 0) {
        index += n;
    }
    if (index < 0 || index >= n) {
        throw std::out_of_range("Frame index " + std::to_string(index) + " is out of range for a trajectory of " +
                                std::to_string(n) + " frames");
    }
    m_current = static_cast<std::size_t>(index);
    return m_frames[m_current];
}

const Frame& MemoryTrajectory::current() const {
    if (m_frames.empty()) {
        throw std::out_of_range("the trajectory has no frames");
    }
    return m_frames[m_current];
}

}  // namespace mdbridge::core
