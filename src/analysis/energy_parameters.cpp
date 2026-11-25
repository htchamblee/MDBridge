// Lennard-Jones parameters and partial charges of the energy analysis without parm7 parameters
// (docs/DIFFERENCES.md AN-25 - AN-27, CIT-04). MDATK v0.0.4.1 used one ad-hoc value per element.
//
// Per residue of the analysed atoms:
//  1. Residue templates of Amber's force fields as OpenMM converts them (ff14SB, OL15, OL3, TIP3P
//     water and ions; src/analysis/energy_lj_tables.inc, tools/gen_energy_lj.py). A residue whose
//     every atom maps to one template atom (by name, then by pdbNames.xml's aliases) takes the
//     template's Lennard-Jones parameters; when the atoms map onto the whole template (hydrogens
//     included) also its partial charges. Amber's terminal variants (N/C amino acids, 5'/3'/N
//     nucleotides) are told apart by their atoms, or, for a residue with missing atoms, by whether
//     it is bonded to the residues before and after it.
//  2. Every other atom: one Amber atom type (parm10.dat, Joung-Cheatham and Li-Merz TIP3P ions)
//     by its element, hydrogens by the element of the atom they are bonded to, carbon and oxygen
//     by their neighbours where bonds are known (the "typed" layer; the design of MDATK's later
//     generic Lennard-Jones tables, same owner).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::analysis {

namespace {

struct AmberType {
    std::string_view name;
    double rstar;
    double epsilon;
    std::string_view reference;
};
struct TemplateAtom {
    std::string_view name;
    std::string_view element;
    std::string_view amberClass;
    double charge;
    double rstar;
    double epsilon;
};
struct ShortCitation {
    std::string_view key;
    std::string_view text;
};
struct Template {
    std::string_view name;
    std::string_view forceField;
    std::string_view family;  // protein, nucleic, water, ion
    std::string_view base;
    std::string_view variant;  // "", N / C (protein), 5 / 3 / N (nucleic)
    int firstAtom;
    int atomCount;
    std::string_view genericAliases;
    std::string_view specificAliases;
    unsigned referenceMask;
};
struct NamePair {
    std::string_view name;
    std::string_view target;
};
struct AtomAlias {
    std::string_view group;
    std::string_view alias;
    std::string_view atom;
};

#include "energy_lj_tables.inc"

constexpr double kHydrogenReach = 1.45;  // Å: a hydrogen without topology bonds belongs to the
                                         // nearest heavy atom of its residue within this distance

const AmberType& amberType(std::string_view name) {
    for (const AmberType& t : kAmberTypes) {
        if (t.name == name) {
            return t;
        }
    }
    throw std::logic_error("energy_lj_tables.inc has no Amber type " + std::string(name));
}

// typed layer: element -> Amber type of an atom of a molecule (parm10.dat) and of a monatomic ion
constexpr std::array<std::array<std::string_view, 2>, 9> kMoleculeTypes{{
    {"C", "C"}, {"N", "N"}, {"O", "O"}, {"P", "P"}, {"S", "S"}, {"F", "F"}, {"CL", "Cl"}, {"BR", "Br"}, {"I", "I"}}};
constexpr std::array<std::array<std::string_view, 2>, 20> kIonTypes{{
    {"LI", "Li+"},  {"NA", "Na+"},  {"K", "K+"},    {"RB", "Rb+"},  {"CS", "Cs+"},  {"F", "F-"},   {"CL", "Cl-"},
    {"BR", "Br-"},  {"I", "I-"},    {"MG", "Mg2+"}, {"CA", "Ca2+"}, {"SR", "Sr2+"}, {"BA", "Ba2+"}, {"MN", "Mn2+"},
    {"FE", "Fe2+"}, {"CO", "Co2+"}, {"NI", "Ni2+"}, {"CU", "Cu2+"}, {"ZN", "Zn2+"}, {"CD", "Cd2+"}}};

template <std::size_t N>
std::optional<std::string_view> lookup(const std::array<std::array<std::string_view, 2>, N>& table,
                                       std::string_view element) {
    for (const auto& row : table) {
        if (row[0] == element) {
            return row[1];
        }
    }
    return std::nullopt;
}

std::string upperAscii(std::string text) {
    for (char& c : text) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return text;
}

// The name without trailing charge signs and digits ("Na+" -> "Na", "Mg2+" -> "Mg").
std::string stripCharge(std::string text) {
    while (!text.empty() && (text.back() == '+' || text.back() == '-' || (text.back() >= '0' && text.back() <= '9'))) {
        text.pop_back();
    }
    return text;
}

std::optional<std::size_t> findTemplate(std::string_view name) {
    for (std::size_t k = 0; k < std::size(kTemplates); ++k) {
        if (kTemplates[k].name == name) {
            return k;
        }
    }
    return std::nullopt;
}

std::span<const TemplateAtom> templateAtoms(const Template& t) {
    return std::span<const TemplateAtom>(kTemplateAtoms).subspan(static_cast<std::size_t>(t.firstAtom),
                                                                 static_cast<std::size_t>(t.atomCount));
}

std::optional<std::string_view> aliasTarget(std::string_view group, std::string_view alias) {
    if (group.empty()) {
        return std::nullopt;
    }
    for (const AtomAlias& a : kAtomAliases) {
        if (a.group == group && a.alias == alias) {
            return a.atom;
        }
    }
    return std::nullopt;
}

bool hasVariants(const Template& t) {
    for (const Template& other : kTemplates) {
        if (!other.variant.empty() && other.base == t.name && other.forceField == t.forceField) {
            return true;
        }
    }
    return false;
}

// The template and, when it is the middle form of a residue with terminal variants, those variants
// (file order).
void addWithVariants(std::size_t k, std::vector<std::size_t>& out, bool& filterable) {
    out.push_back(k);
    const Template& t = kTemplates[k];
    if (!t.variant.empty() || !hasVariants(t)) {
        return;
    }
    filterable = true;
    for (std::size_t v = 0; v < std::size(kTemplates); ++v) {
        if (!kTemplates[v].variant.empty() && kTemplates[v].base == t.name && kTemplates[v].forceField == t.forceField) {
            out.push_back(v);
        }
    }
}

// The template atom of each residue atom (-1: none), or std::nullopt when some atom has none or
// two atoms take one template atom. Names first, then the template's residue-specific and
// generic aliases (pdbNames.xml); an atom's known element must be the template atom's.
std::optional<std::vector<int>> mapAtoms(const Template& t, std::span<const std::string> names,
                                         std::span<const std::string> elements, bool& complete) {
    const std::span<const TemplateAtom> atoms = templateAtoms(t);
    std::vector<int> map(names.size(), -1);
    std::vector<bool> used(atoms.size(), false);
    const auto indexOf = [&](std::string_view name) -> int {
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            if (atoms[a].name == name) {
                return static_cast<int>(a);
            }
        }
        return -1;
    };
    for (std::size_t k = 0; k < names.size(); ++k) {
        const int a = indexOf(names[k]);
        if (a >= 0) {
            if (used[static_cast<std::size_t>(a)]) {
                return std::nullopt;  // a repeated atom name
            }
            used[static_cast<std::size_t>(a)] = true;
            map[k] = a;
        }
    }
    for (std::size_t k = 0; k < names.size(); ++k) {
        if (map[k] >= 0) {
            continue;
        }
        std::optional<std::string_view> target = aliasTarget(t.specificAliases, names[k]);
        if (!target) {
            target = aliasTarget(t.genericAliases, names[k]);
        }
        const int a = target ? indexOf(*target) : -1;
        if (a < 0 || used[static_cast<std::size_t>(a)]) {
            return std::nullopt;
        }
        used[static_cast<std::size_t>(a)] = true;
        map[k] = a;
    }
    for (std::size_t k = 0; k < names.size(); ++k) {
        const std::string_view element = atoms[static_cast<std::size_t>(map[k])].element;
        if (!elements[k].empty() && !element.empty() && elements[k] != element) {
            return std::nullopt;
        }
    }
    complete = std::all_of(used.begin(), used.end(), [](bool u) { return u; });
    return map;
}

class Assigner {
public:
    Assigner(const core::AtomGroup& atoms, const core::Frame& frame)
        : m_atoms(atoms), m_t(atoms.topology()), m_frame(frame) {
        m_residueHasTopologyBonds.assign(m_t.residueCount(), false);
        m_topologyPartners.resize(m_t.atomCount());
        if (m_t.bonds) {
            for (const core::Bond& b : *m_t.bonds) {
                m_topologyPartners[static_cast<std::size_t>(b[0])].push_back(b[1]);
                m_topologyPartners[static_cast<std::size_t>(b[1])].push_back(b[0]);
                m_residueHasTopologyBonds[static_cast<std::size_t>(m_t.atomResidue[static_cast<std::size_t>(b[0])])] = true;
                m_residueHasTopologyBonds[static_cast<std::size_t>(m_t.atomResidue[static_cast<std::size_t>(b[1])])] = true;
            }
        }
        m_elements = core::atomElements(core::AtomGroup(atoms.sharedTopology()));
        m_extraPoint.assign(m_t.atomCount(), false);
        for (std::size_t a = 0; a < m_t.atomCount(); ++a) {
            m_extraPoint[a] = isExtraPoint(static_cast<std::int32_t>(a));
        }
    }

    EnergyParameters run() {
        EnergyParameters out;
        const std::size_t n = m_atoms.size();
        out.rstar.assign(n, 0.0);
        out.epsilon.assign(n, 0.0);
        out.charge.assign(n, 0.0);
        out.layer.assign(n, EnergyParameters::Layer::Typed);
        out.templateCharge.assign(n, false);
        // per atom of the topology: template index and template atom, the template's charges usable
        std::vector<int> atomTemplate(m_t.atomCount(), -1);
        std::vector<int> atomTemplateAtom(m_t.atomCount(), -1);
        std::vector<bool> atomCharged(m_t.atomCount(), false);
        std::vector<bool> residueDone(m_t.residueCount(), false);
        for (const std::int32_t atom : m_atoms.indices()) {
            const std::int32_t r = m_t.atomResidue[static_cast<std::size_t>(atom)];
            if (residueDone[static_cast<std::size_t>(r)]) {
                continue;
            }
            residueDone[static_cast<std::size_t>(r)] = true;
            matchResidue(r, atomTemplate, atomTemplateAtom, atomCharged);
        }
        std::vector<bool> chargelessResidue(m_t.residueCount(), false);
        for (std::size_t k = 0; k < n; ++k) {
            const std::int32_t atom = m_atoms[k];
            const auto ua = static_cast<std::size_t>(atom);
            if (atomTemplate[ua] >= 0) {
                const Template& t = kTemplates[static_cast<std::size_t>(atomTemplate[ua])];
                const TemplateAtom& ta = templateAtoms(t)[static_cast<std::size_t>(atomTemplateAtom[ua])];
                out.rstar[k] = ta.rstar;
                out.epsilon[k] = ta.epsilon;
                out.layer[k] = EnergyParameters::Layer::Template;
                if (atomCharged[ua]) {
                    out.charge[k] = ta.charge;
                    out.templateCharge[k] = true;
                }
                addUnique(out.forceFields, std::string(t.forceField));
                for (std::size_t bit = 0; bit < std::size(kReferenceKeys); ++bit) {
                    if ((t.referenceMask >> bit) & 1u) {
                        addUnique(out.references, std::string(kReferenceKeys[bit]));
                    }
                }
            } else {
                const AmberType& type = typed(atom, out);
                out.rstar[k] = type.rstar;
                out.epsilon[k] = type.epsilon;
                addUnique(out.references, std::string(type.reference));
            }
            if (!out.templateCharge[k]) {
                chargelessResidue[static_cast<std::size_t>(m_t.atomResidue[ua])] = true;
            }
        }
        out.chargelessResidues = static_cast<std::size_t>(
            std::count(chargelessResidue.begin(), chargelessResidue.end(), true));
        out.guessedBonds = m_guessed;
        return out;
    }

private:
    static void addUnique(std::vector<std::string>& list, const std::string& item) {
        if (!item.empty() && std::find(list.begin(), list.end(), item) == list.end()) {
            list.push_back(item);
        }
    }

    // AN-25: the template of residue r, if one matches (see the file comment).
    void matchResidue(std::int32_t r, std::vector<int>& atomTemplate, std::vector<int>& atomTemplateAtom,
                      std::vector<bool>& atomCharged) {
        // extra points take no part (AN-26: Amber's EP, no Lennard-Jones term); a residue with one
        // never takes a template's charges (a 4-site water is no TIP3P water)
        std::vector<std::int32_t> atoms;
        bool extraPoints = false;
        for (const std::int32_t a : m_t.residueAtoms(r)) {
            if (m_extraPoint[static_cast<std::size_t>(a)]) {
                extraPoints = true;
            } else {
                atoms.push_back(a);
            }
        }
        if (atoms.empty()) {
            return;
        }
        std::vector<std::string> names;
        std::vector<std::string> elements;
        for (const std::int32_t a : atoms) {
            names.push_back(m_t.name(a));
            elements.push_back(m_elements[static_cast<std::size_t>(a)]);
        }
        const auto assign = [&](std::size_t k, const std::vector<int>& map, bool charged) {
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const auto ua = static_cast<std::size_t>(atoms[i]);
                atomTemplate[ua] = static_cast<int>(k);
                atomTemplateAtom[ua] = map[i];
                atomCharged[ua] = charged && !extraPoints;
            }
        };
        const std::string resname = upperAscii(core::pythonStrip(m_t.resname(r)));
        if (atoms.size() == 1 && !extraPoints) {
            if (const std::optional<std::size_t> ion = matchIon(r, elements[0])) {
                assign(*ion, {0}, true);
                return;
            }
        }
        std::vector<std::size_t> candidates;
        bool filterable = false;
        if (const std::optional<std::size_t> k = findTemplate(resname); k && kTemplates[*k].family != "ion") {
            addWithVariants(*k, candidates, filterable);
        } else {
            for (const NamePair& alias : kResidueAliases) {
                if (alias.name == resname) {
                    if (const std::optional<std::size_t> target = findTemplate(alias.target)) {
                        addWithVariants(*target, candidates, filterable);
                    }
                }
            }
        }
        if (candidates.empty()) {
            return;
        }
        struct Match {
            std::size_t k;
            std::vector<int> map;
            bool complete;
        };
        std::vector<Match> containing;
        for (const std::size_t k : candidates) {
            bool complete = false;
            if (std::optional<std::vector<int>> map = mapAtoms(kTemplates[k], names, elements, complete)) {
                containing.push_back({k, std::move(*map), complete});
            }
        }
        const auto completes = std::count_if(containing.begin(), containing.end(), [](const Match& m) { return m.complete; });
        if (completes == 1) {
            for (const Match& m : containing) {
                if (m.complete) {
                    assign(m.k, m.map, true);
                }
            }
            return;
        }
        if (filterable) {
            const bool start = !linked(r - 1, r);
            const bool end = !linked(r, r + 1);
            std::erase_if(containing, [&](const Match& m) {
                const Template& t = kTemplates[m.k];
                if (t.variant.empty() && !hasVariants(t)) {
                    return false;  // no terminal forms (caps): any position
                }
                std::string_view wanted;
                if (t.family == "protein") {
                    wanted = start && end ? "-" : start ? "N" : end ? "C" : "";
                } else {
                    wanted = start && end ? "N" : start ? "5" : end ? "3" : "";
                }
                return t.variant != wanted;
            });
        }
        if (containing.empty()) {
            return;
        }
        // the Lennard-Jones parameters of every atom must agree between the remaining templates
        for (const Match& m : containing) {
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const TemplateAtom& a = templateAtoms(kTemplates[m.k])[static_cast<std::size_t>(m.map[i])];
                const TemplateAtom& b =
                    templateAtoms(kTemplates[containing[0].k])[static_cast<std::size_t>(containing[0].map[i])];
                if (a.rstar != b.rstar || a.epsilon != b.epsilon) {
                    return;
                }
            }
        }
        assign(containing[0].k, containing[0].map, false);
    }

    // A one-atom residue: the ion template named like the residue (as written, upper-cased,
    // without charge signs and digits) or else like the atom, whose element is the atom's.
    std::optional<std::size_t> matchIon(std::int32_t r, const std::string& element) const {
        const std::int32_t atom = m_t.residueAtoms(r)[0];
        const std::string resname = core::pythonStrip(m_t.resname(r));
        const std::string name = core::pythonStrip(m_t.name(atom));
        for (const std::string& key : {resname, upperAscii(resname), upperAscii(stripCharge(resname)), name,
                                       upperAscii(name), upperAscii(stripCharge(name))}) {
            if (key.empty()) {
                continue;
            }
            if (const std::optional<std::size_t> k = findTemplate(key)) {
                const Template& t = kTemplates[*k];
                if (t.family == "ion" && templateAtoms(t)[0].element == element) {
                    return k;
                }
            }
        }
        return std::nullopt;
    }

    // Residues r1 and r2 (= r1 + 1) continue one chain. A chain breaks between them only when the
    // evidence agrees: a new chain (another segment, or another chain ID where the topology has
    // chain IDs), or a gap or reset in the residue numbering (resid(r2) is neither resid(r1) + 1 nor
    // resid(r1) with another insertion code) without a bond between them (the topology's bonds
    // when either residue has any, else the bonds guessed on the frame, AN-14). So consecutive
    // numbering within one chain links residues whatever their geometry, and a bond bridges a
    // numbering gap.
    bool linked(std::int32_t r1, std::int32_t r2) {
        if (r1 < 0 || r2 < 0 || static_cast<std::size_t>(r2) >= m_t.residueCount() ||
            static_cast<std::size_t>(r1) >= m_t.residueCount()) {
            return false;
        }
        const auto u1 = static_cast<std::size_t>(r1);
        const auto u2 = static_cast<std::size_t>(r2);
        if (m_t.residueSegment[u1] != m_t.residueSegment[u2]) {
            return false;
        }
        if (m_t.chainIDs) {
            const auto first = [&](std::int32_t r) { return static_cast<std::size_t>(m_t.residueAtoms(r)[0]); };
            if ((*m_t.chainIDs)[first(r1)] != (*m_t.chainIDs)[first(r2)]) {
                return false;
            }
        }
        const std::int64_t n1 = m_t.resids[u1];
        const std::int64_t n2 = m_t.resids[u2];
        const bool otherInsertion = m_t.icodes && (*m_t.icodes)[u1] != (*m_t.icodes)[u2];
        if (n2 == n1 + 1 || (n2 == n1 && otherInsertion)) {
            return true;
        }
        const bool topology = m_residueHasTopologyBonds[u1] || m_residueHasTopologyBonds[u2];
        for (const std::int32_t a : m_t.residueAtoms(r1)) {
            for (const std::int32_t b : topology ? m_topologyPartners[static_cast<std::size_t>(a)] : guessedPartners(a)) {
                if (m_t.atomResidue[static_cast<std::size_t>(b)] == r2) {
                    return true;
                }
            }
        }
        return false;
    }

    // Bonds guessed on the frame among the atoms of the analysed residues and their neighbours
    // (computed once, on first use).
    const std::vector<std::int32_t>& guessedPartners(std::int32_t atom) {
        if (!m_guessed) {
            m_guessed = true;
            m_guessedPartners.resize(m_t.atomCount());
            std::vector<bool> keep(m_t.residueCount(), false);
            for (const std::int32_t a : m_atoms.indices()) {
                const std::int32_t r = m_t.atomResidue[static_cast<std::size_t>(a)];
                for (std::int32_t q = r - 1; q <= r + 1; ++q) {
                    if (q >= 0 && static_cast<std::size_t>(q) < m_t.residueCount()) {
                        keep[static_cast<std::size_t>(q)] = true;
                    }
                }
            }
            std::vector<std::int32_t> group;
            for (std::size_t a = 0; a < m_t.atomCount(); ++a) {
                if (keep[static_cast<std::size_t>(m_t.atomResidue[a])]) {
                    group.push_back(static_cast<std::int32_t>(a));
                }
            }
            const core::AtomGroup atoms(m_atoms.sharedTopology(), group);
            const std::vector<std::int64_t> pairs = guessBonds(atoms, m_frame);
            for (std::size_t p = 0; p + 1 < pairs.size(); p += 2) {
                const std::int32_t i = group[static_cast<std::size_t>(pairs[p])];
                const std::int32_t j = group[static_cast<std::size_t>(pairs[p + 1])];
                m_guessedPartners[static_cast<std::size_t>(i)].push_back(j);
                m_guessedPartners[static_cast<std::size_t>(j)].push_back(i);
            }
        }
        return m_guessedPartners[static_cast<std::size_t>(atom)];
    }

    // The bonded neighbours of a heavy atom: the topology's when its residue has topology bonds,
    // else the guessed ones; extra points are no neighbours.
    std::vector<std::int32_t> neighbours(std::int32_t atom) {
        const std::int32_t r = m_t.atomResidue[static_cast<std::size_t>(atom)];
        std::vector<std::int32_t> out = m_residueHasTopologyBonds[static_cast<std::size_t>(r)]
                                            ? m_topologyPartners[static_cast<std::size_t>(atom)]
                                            : guessedPartners(atom);
        std::erase_if(out, [&](std::int32_t b) { return m_extraPoint[static_cast<std::size_t>(b)]; });
        return out;
    }

    // An extra point (virtual site, lone pair): an atom whose element is no chemical element and
    // whose name is a water virtual site of OpenMM's force fields (M, M1, M2, LP1, LP2, OM), MW
    // (GROMACS' TIP4P) or starts with EP or LP (Amber's EPW, EP1, ...; lone pairs); or an atom named
    // like the virtual site of its residue in those force fields (OM of TP4E, whose name reads as
    // oxygen).
    bool isExtraPoint(std::int32_t atom) const {
        const std::string name = upperAscii(core::pythonStrip(m_t.name(atom)));
        const std::string resname = upperAscii(core::pythonStrip(m_t.resname(m_t.atomResidue[static_cast<std::size_t>(atom)])));
        const auto pair = [](const NamePair& p, const NamePair& q) {
            return p.name != q.name ? p.name < q.name : p.target < q.target;
        };
        if (std::binary_search(std::begin(kVirtualSites), std::end(kVirtualSites), NamePair{resname, name}, pair)) {
            return true;
        }
        if (std::binary_search(std::begin(kElementSymbols), std::end(kElementSymbols),
                               std::string_view(m_elements[static_cast<std::size_t>(atom)]))) {
            return false;
        }
        if (name == "MW" || name.rfind("EP", 0) == 0 || name.rfind("LP", 0) == 0) {
            return true;
        }
        return std::any_of(std::begin(kVirtualSites), std::end(kVirtualSites),
                           [&](const NamePair& v) { return v.target == name; });
    }

    // The atom a hydrogen is bonded to: by the topology's bonds (its first partner), else the
    // nearest non-hydrogen atom of its residue within 1.45 Å on the frame; -1 for none.
    std::int32_t hydrogenParent(std::int32_t h) const {
        std::int32_t bonded = -1;
        for (const std::int32_t b : m_topologyPartners[static_cast<std::size_t>(h)]) {
            if (!m_extraPoint[static_cast<std::size_t>(b)] && (bonded < 0 || b < bonded)) {
                bonded = b;
            }
        }
        if (bonded >= 0) {
            return bonded;
        }
        const std::int32_t r = m_t.atomResidue[static_cast<std::size_t>(h)];
        const auto at = [&](std::int32_t a, std::size_t c) {
            return static_cast<double>(m_frame.positions[static_cast<std::size_t>(a) * 3 + c]);
        };
        std::int32_t best = -1;
        double bestDistance = 0.0;
        for (const std::int32_t a : m_t.residueAtoms(r)) {
            if (a == h || m_elements[static_cast<std::size_t>(a)] == "H" || m_extraPoint[static_cast<std::size_t>(a)]) {
                continue;
            }
            const double dx = at(a, 0) - at(h, 0);
            const double dy = at(a, 1) - at(h, 1);
            const double dz = at(a, 2) - at(h, 2);
            const double d = std::sqrt((dx * dx + dy * dy) + dz * dz);
            if (d <= kHydrogenReach && (best < 0 || d < bestDistance)) {
                best = a;
                bestDistance = d;
            }
        }
        return best;
    }

    // AN-26: the Amber type of an atom no template covers.
    const AmberType& typed(std::int32_t atom, EnergyParameters& out) {
        const std::string& element = m_elements[static_cast<std::size_t>(atom)];
        if (m_extraPoint[static_cast<std::size_t>(atom)]) {
            return amberType("EP");  // no Lennard-Jones term
        }
        if (element == "H") {
            const std::int32_t parent = hydrogenParent(atom);
            const std::string parentElement = parent >= 0 ? m_elements[static_cast<std::size_t>(parent)] : "N";
            return amberType(parentElement == "N" ? "H" : parentElement == "O" ? "HO" : parentElement == "S" ? "HS" : "H1");
        }
        const std::int32_t r = m_t.atomResidue[static_cast<std::size_t>(atom)];
        if (m_t.residueAtoms(r).size() == 1) {
            if (const std::optional<std::string_view> ion = lookup(kIonTypes, element)) {
                return amberType(*ion);
            }
        }
        if (const std::optional<std::string_view> type = lookup(kMoleculeTypes, element)) {
            const std::vector<std::int32_t> bonded =
                element == "C" || element == "O" ? neighbours(atom) : std::vector<std::int32_t>{};
            if (element == "C" && bonded.size() == 4) {
                return amberType("CT");  // sp3 carbon
            }
            if (element == "O" && bonded.size() == 2) {
                std::size_t hydrogens = 0;
                for (const std::int32_t b : bonded) {
                    if (m_elements[static_cast<std::size_t>(b)] == "H") {
                        ++hydrogens;
                    }
                }
                return amberType(hydrogens == 2 ? "OW" : hydrogens == 1 ? "OH" : "OS");
            }
            return amberType(*type);
        }
        if (const std::optional<std::string_view> ion = lookup(kIonTypes, element)) {
            return amberType(*ion);  // a metal inside a larger residue: its ion's parameters
        }
        if (std::find(out.defaultElements.begin(), out.defaultElements.end(), element) == out.defaultElements.end()) {
            out.defaultElements.push_back(element);
        }
        ++out.defaultAtoms;
        return amberType("C");
    }

    const core::AtomGroup& m_atoms;
    const core::Topology& m_t;
    const core::Frame& m_frame;
    std::vector<std::string> m_elements;  // every atom of the topology
    std::vector<bool> m_extraPoint;       // every atom of the topology
    std::vector<bool> m_residueHasTopologyBonds;
    std::vector<std::vector<std::int32_t>> m_topologyPartners;
    bool m_guessed = false;
    std::vector<std::vector<std::int32_t>> m_guessedPartners;
};

}  // namespace

EnergyParameters energyParameters(const core::AtomGroup& atoms, const core::Frame& frame) {
    if (atoms.empty()) {
        return {};
    }
    return Assigner(atoms, frame).run();
}

std::string energyShortCitation(std::string_view key) {
    for (const ShortCitation& c : kShortCitations) {
        if (c.key == key) {
            return std::string(c.text);
        }
    }
    return std::string(key);
}

}  // namespace mdbridge::analysis
