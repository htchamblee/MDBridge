// Restraint dictionaries of MDATK's mdatk/dictionaries.py (MDATK's own code, ported): COMPONENTS,
// smiles_mapping / check_smiles / identify_smiles and finish_dictionary / rename_dictionary.

#include "mdbridge/chem/dictionaries.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <set>
#include <stdexcept>

#include "mdbridge/chem/reference_tables.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::chem {

namespace {

[[noreturn]] void fail(const std::string& message) { throw DictionaryError(message); }

std::string join(const std::vector<std::string>& items, std::string_view sep) {
    std::string out;
    for (std::size_t k = 0; k < items.size(); ++k) {
        if (k != 0) {
            out += sep;
        }
        out += items[k];
    }
    return out;
}

bool contains(const std::vector<std::string>& items, std::string_view value) {
    return std::find(items.begin(), items.end(), value) != items.end();
}

std::string repr(const std::string& text) { return core::pythonStringRepr(text); }

// str.upper() of a residue code, ASCII letters only (DIFFERENCES DIC-02).
std::string asciiUpper(std::string text) {
    for (char& c : text) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 32);
        }
    }
    return text;
}

bool startsWith(std::string_view text, std::string_view prefix) {
    return text.substr(0, prefix.size()) == prefix;
}

// --- reference data -----------------------------------------------------------------------------

// MDATK's own data of the four parents (dictionaries.py _PARENT_DATA): the nucleoside names and
// the isomeric L SMILES (deprotonated phosphate, Kekule form); atoms, bonds, hydrogens and naming
// rules come from the CCD (reference_tables.hpp).
struct ParentData {
    std::string_view code;
    std::string_view nucleoside;
    std::string_view smiles;
    std::string_view smiles5;
};

constexpr ParentData kParents[] = {
    {"0A", "L-ADENOSINE", "[O-]P(=O)([O-])OC[C@@H]1O[C@H](N2C=NC3=C(N)N=CN=C32)[C@@H](O)[C@H]1O",
     "OC[C@@H]1O[C@H](N2C=NC3=C(N)N=CN=C32)[C@@H](O)[C@H]1O"},
    {"0C", "L-CYTIDINE", "[O-]P(=O)([O-])OC[C@@H]1O[C@H](N2C=CC(N)=NC2=O)[C@@H](O)[C@H]1O",
     "OC[C@@H]1O[C@H](N2C=CC(N)=NC2=O)[C@@H](O)[C@H]1O"},
    {"0G", "L-GUANOSINE", "[O-]P(=O)([O-])OC[C@@H]1O[C@H](N2C=NC3=C2N=C(N)NC3=O)[C@@H](O)[C@H]1O",
     "OC[C@@H]1O[C@H](N2C=NC3=C2N=C(N)NC3=O)[C@@H](O)[C@H]1O"},
    {"0U", "L-URIDINE", "[O-]P(=O)([O-])OC[C@@H]1O[C@H](N2C=CC(=O)NC2=O)[C@@H](O)[C@H]1O",
     "OC[C@@H]1O[C@H](N2C=CC(=O)NC2=O)[C@@H](O)[C@H]1O"},
};

constexpr std::string_view kPhosphate[] = {"OP3", "P", "OP1", "OP2"};

bool isPhosphate(std::string_view name) {
    return std::find(std::begin(kPhosphate), std::end(kPhosphate), name) != std::end(kPhosphate);
}

// PROTONATION_SITES: heavy atoms whose hydrogen count may differ (phosphate) -> the H name
std::optional<std::string> protonationName(std::string_view heavy) {
    if (heavy == "OP1" || heavy == "OP2" || heavy == "OP3") {
        return "H" + std::string(heavy);
    }
    return std::nullopt;
}

std::vector<Component> buildComponents() {
    std::map<std::string, Component> comps;
    std::vector<std::string> parents;
    for (const reference::ComponentInfo& info : reference::components()) {
        const std::string code(info.code);
        const auto data = std::find_if(std::begin(kParents), std::end(kParents),
                                       [&](const ParentData& p) { return p.code == info.code; });
        if (data == std::end(kParents)) {
            throw std::logic_error("no SMILES for component " + code);
        }
        Component c;
        c.code = code;
        c.parent = code;
        c.base = std::string(1, info.base);
        c.name = std::string(info.name);
        c.smiles = std::string(data->smiles);
        for (const auto& a : reference::componentAtoms(info.code)) {
            c.atoms.push_back({std::string(a.name), std::string(a.element), a.charge});
        }
        for (const auto& b : reference::componentBonds(info.code)) {
            c.bonds.push_back({std::string(b.atom1), std::string(b.atom2), b.order});
        }
        for (const auto& h : reference::componentHydrogens(info.code)) {
            if (c.hydrogens.empty() || c.hydrogens.back().first != h.heavy) {
                c.hydrogens.emplace_back(std::string(h.heavy), std::vector<std::string>{});
            }
            c.hydrogens.back().second.emplace_back(h.hydrogen);
        }
        for (const auto& r : reference::hydrogenPositions(info.code)) {
            c.hpos.push_back({std::string(r.hydrogen), std::string(r.x), std::string(r.a), std::string(r.b),
                              std::string(r.kind)});
        }
        c.description = "CCD entry (internal residue; OP3 and HO3' are leaving atoms of the polymer link)";

        Component five = c;
        five.code = code + "5";
        five.terminal = "5";
        five.name = std::string(data->nucleoside) + " (5'-TERMINAL, 5'-OH)";
        five.smiles = std::string(data->smiles5);
        std::erase_if(five.atoms, [](const ComponentAtomName& a) { return isPhosphate(a.name); });
        std::erase_if(five.bonds,
                      [](const ComponentBondRow& b) { return isPhosphate(b.atom1) || isPhosphate(b.atom2); });
        five.hydrogens.insert(five.hydrogens.begin(), {"O5'", {"HO5'"}});
        five.description = "5'-terminal variant: no 5'-phosphate, HO5' on O5'";

        Component three = c;
        three.code = code + "3";
        three.terminal = "3";
        three.name = c.name + " (3'-TERMINAL, NO OP3)";
        three.leaving = {"OP3"};
        three.description = "3'-terminal variant: P, OP1, OP2 without OP3, free 3'-OH (generated from the parent; "
                            "OP3 is removed after renaming)";

        parents.push_back(code);
        comps[five.code] = std::move(five);
        comps[three.code] = std::move(three);
        comps[code] = std::move(c);
    }
    std::vector<Component> out;
    for (const std::string& p : parents) {
        out.push_back(comps.at(p + "5"));
        out.push_back(comps.at(p + "3"));
    }
    for (const std::string& p : parents) {
        out.push_back(comps.at(p));
    }
    return out;
}

// --- geometry (pure Python arithmetic in MDATK's order) ---------------------------------------

using Vec = std::array<double, 3>;

Vec sub(const Vec& a, const Vec& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec cross(const Vec& a, const Vec& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Vec& a, const Vec& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

double torsion(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3) {
    const Vec b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
    double n1 = std::sqrt(dot(b1, b1));
    if (n1 == 0.0) {
        n1 = 1.0;
    }
    const Vec b1n = {b1[0] / n1, b1[1] / n1, b1[2] / n1};
    const double d0 = dot(b0, b1n), d2 = dot(b2, b1n);
    const Vec v = sub(b0, Vec{d0 * b1n[0], d0 * b1n[1], d0 * b1n[2]});
    const Vec w = sub(b2, Vec{d2 * b1n[0], d2 * b1n[1], d2 * b1n[2]});
    constexpr double kDegrees = 180.0 / 3.141592653589793;
    return std::atan2(dot(cross(b1n, v), w), dot(v, w)) * kDegrees;
}

double volume(const Vec& c, const Vec& a, const Vec& b, const Vec& d) {
    const Vec va = sub(a, c), vb = sub(b, c), vd = sub(d, c);
    return va[0] * (vb[1] * vd[2] - vb[2] * vd[1]) - va[1] * (vb[0] * vd[2] - vb[2] * vd[0]) +
           va[2] * (vb[0] * vd[1] - vb[1] * vd[0]);
}

// _parity: +1 if seq is an even permutation of ref, -1 if odd (every element of seq is in ref).
int parity(const std::vector<std::string>& seq, const std::vector<std::string>& ref) {
    std::vector<std::ptrdiff_t> p;
    for (const std::string& x : seq) {
        p.push_back(std::find(ref.begin(), ref.end(), x) - ref.begin());
    }
    int inversions = 0;
    for (std::size_t a = 0; a < p.size(); ++a) {
        for (std::size_t b = a + 1; b < p.size(); ++b) {
            if (p[a] > p[b]) {
                ++inversions;
            }
        }
    }
    return inversions % 2 != 0 ? -1 : 1;
}

int expectedSign(const std::vector<std::string>& refOrder, int refSign, const std::vector<std::string>& order) {
    return refSign * parity(order, refOrder);
}

std::vector<std::string> orderNames(const ChiralityReference& c) {
    return {c.order.begin(), c.order.end()};
}

// An insertion-ordered str -> value dict (Python's dict: setting a key again keeps its place).
template <typename T>
class OrderedMap {
public:
    void set(const std::string& key, T value) {
        const auto it = m_index.find(key);
        if (it != m_index.end()) {
            m_items[it->second].second = std::move(value);
            return;
        }
        m_index[key] = m_items.size();
        m_items.emplace_back(key, std::move(value));
    }
    const T* get(const std::string& key) const {
        const auto it = m_index.find(key);
        return it == m_index.end() ? nullptr : &m_items[it->second].second;
    }
    bool has(const std::string& key) const { return m_index.count(key) != 0; }
    const std::vector<std::pair<std::string, T>>& items() const { return m_items; }
    std::size_t size() const { return m_items.size(); }
    bool empty() const { return m_items.empty(); }

private:
    std::vector<std::pair<std::string, T>> m_items;
    std::map<std::string, std::size_t> m_index;
};

using Coords = OrderedMap<Vec>;
using NameMap = OrderedMap<std::string>;

std::string mapped(const NameMap& mapping, const std::string& name) {
    const std::string* v = mapping.get(name);
    return v ? *v : name;
}

}  // namespace

// --- components -----------------------------------------------------------------------------------

std::span<const std::string_view> dictionaryGroups() {
    static constexpr std::string_view kGroups[] = {"RNA", "non-polymer"};
    return kGroups;
}

const std::vector<std::string>& Component::hydrogenNames(std::string_view heavy) const {
    static const std::vector<std::string> kNone;
    const std::vector<std::string>* found = &kNone;
    for (const auto& [x, hs] : hydrogens) {  // dict(self.hydrogens): the last entry of a name wins
        if (x == heavy) {
            found = &hs;
        }
    }
    return *found;
}

bool Component::isLeaving(std::string_view atom) const { return contains(leaving, atom); }

std::vector<std::string> Component::residueAtoms() const {
    std::vector<std::string> out;
    for (const auto& a : atoms) {
        if (!isLeaving(a.name)) {
            out.push_back(a.name);
        }
    }
    for (const auto& [x, hs] : hydrogens) {
        if (!isLeaving(x)) {
            out.insert(out.end(), hs.begin(), hs.end());
        }
    }
    return out;
}

std::string Component::glycosidicN() const { return (base == "A" || base == "G") ? "N9" : "N1"; }

std::vector<ChiralityReference> Component::chirality() const {
    std::vector<ChiralityReference> out;
    for (const reference::ChiralCentre& row : reference::lRiboseChirality()) {
        ChiralityReference c;
        c.centre = std::string(row.centre);
        const std::string_view nbrs[] = {row.a, row.b, row.c};
        for (std::size_t k = 0; k < 3; ++k) {
            c.order[k] = nbrs[k] == "N" ? glycosidicN() : std::string(nbrs[k]);
        }
        c.order[3] = hydrogenNames(c.centre).at(0);
        c.sign = row.sign;
        out.push_back(std::move(c));
    }
    return out;
}

MolGraph Component::graph(bool withLeaving) const {
    MolGraph g;
    std::map<std::string, int> index;
    for (const auto& a : atoms) {
        if (withLeaving || !isLeaving(a.name)) {
            const std::vector<std::string>& hs = hydrogenNames(a.name);
            index[a.name] = g.addAtom(a.element, a.charge, static_cast<int>(hs.size()), a.name, hs);
        }
    }
    for (const auto& b : bonds) {
        const auto i = index.find(b.atom1), j = index.find(b.atom2);
        if (i != index.end() && j != index.end()) {
            g.addBond(i->second, j->second, b.order);
        }
    }
    return g;
}

std::string Component::formula() const { return graph().formula(); }

const std::vector<Component>& dictionaryComponents() {
    static const std::vector<Component> comps = buildComponents();
    return comps;
}

const Component* findComponent(std::string_view code) {
    for (const Component& c : dictionaryComponents()) {
        if (c.code == code) {
            return &c;
        }
    }
    return nullptr;
}

const Component& getComponent(std::string_view code) {
    if (const Component* c = findComponent(asciiUpper(core::pythonStrip(code)))) {
        return *c;
    }
    std::vector<std::string> known;
    for (const Component& c : dictionaryComponents()) {
        known.push_back(c.code);
    }
    fail("unknown residue code " + repr(std::string(code)) + " (known: " + join(known, ", ") + ")");
}

// --- SMILES against a component ------------------------------------------------------------------

std::string describeMismatch(const Component& comp, const MolGraph& tgt, std::string_view what) {
    const MolGraph ref = comp.graph();
    return std::string(what) + " (" + tgt.formula() + ", " + std::to_string(tgt.size()) +
           " heavy atoms) does not match " + comp.code + " (" + ref.formula() + ", " + std::to_string(ref.size()) +
           " heavy atoms)";
}

SmilesMapping smilesMapping(const Component& comp, std::string_view smiles) {
    const MolGraph tgt = parseSmiles(smiles);
    const MolGraph ref = comp.graph();
    const std::optional<std::vector<int>> m = matchGraph(ref, tgt);
    if (!m) {
        fail(describeMismatch(comp, tgt, "the SMILES"));
    }
    SmilesMapping out;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        out.names[(*m)[i]] = ref.names[i];
    }
    for (std::size_t i = 0; i < ref.size(); ++i) {
        const int got = tgt.hcount[static_cast<std::size_t>((*m)[i])];
        if (ref.hcount[i] != got) {
            if (protonationName(ref.names[i])) {
                out.notes.push_back(ref.names[i] + " has " + std::to_string(got) + " H (reference " +
                                    std::to_string(ref.hcount[i]) + "): different phosphate protonation");
            } else {
                fail("the SMILES has " + std::to_string(got) + " H on " + ref.names[i] + ", " + comp.code + " has " +
                     std::to_string(ref.hcount[i]) + " (different tautomer / protonation)");
            }
        }
    }
    std::map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        index.emplace(ref.names[i], i);  // {name: i}: names are unique
    }
    const std::vector<ChiralityReference> centres = comp.chirality();
    for (const ChiralityReference& c : centres) {
        const int t = (*m)[index.at(c.centre)];
        const auto it = tgt.stereo.find(t);
        if (it == tgt.stereo.end()) {
            fail("the SMILES does not specify the configuration at " + c.centre + " (isomeric SMILES needed)");
        }
        std::vector<std::string> named;
        for (const int x : it->second.order) {
            named.push_back(x == MolGraph::kHydrogen ? "H" : out.names.at(x));
        }
        std::vector<std::string> refNamed = orderNames(c);
        refNamed[3] = "H";
        std::vector<std::string> a = named, b = refNamed;
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        if (a != b) {
            fail("unexpected neighbours of " + c.centre + " in the SMILES");
        }
        const int smilesSign = it->second.tag == "@" ? 1 : -1;  // '@': V(n0, n1, n2) > 0
        if (expectedSign(named, smilesSign, refNamed) != c.sign) {
            fail("the SMILES has the opposite configuration at " + c.centre +
                 " (D sugar or another stereoisomer); " + comp.code + " is an L-nucleotide");
        }
    }
    std::vector<std::string> extra;
    for (const auto& [t, _s] : tgt.stereo) {
        const std::string& name = out.names.at(t);
        if (std::none_of(centres.begin(), centres.end(), [&](const ChiralityReference& c) { return c.centre == name; })) {
            extra.push_back(name);
        }
    }
    if (!extra.empty()) {
        std::sort(extra.begin(), extra.end());
        out.notes.push_back("stereo tags ignored at " + join(extra, ", "));
    }
    return out;
}

std::vector<std::string> checkSmiles(std::string_view code, std::string_view smiles) {
    return smilesMapping(getComponent(code), smiles).notes;
}

std::vector<std::string> identifySmiles(std::string_view smiles) {
    std::vector<std::string> out;
    for (const Component& comp : dictionaryComponents()) {
        try {
            smilesMapping(comp, smiles);
        } catch (const DictionaryError&) {
            continue;
        }
        out.push_back(comp.code);
    }
    return out;
}

// --- finishing a dictionary ------------------------------------------------------------------------

namespace {

using io::CifBlock;
using io::CifDocument;
using io::CifLoop;

// _element: the type symbol, else the first letter of the atom name ("X" without one), as
// "Cl" (first character upper case, the rest lower case).
std::string elementOf(const std::string& symbol, const std::string& name) {
    std::string sym = core::pythonStrip(symbol);
    if (sym.empty() || sym == "." || sym == "?" || sym == ".?") {
        const std::size_t len = name.empty() ? 0 : name.find_first_not_of(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
        sym = (len == 0) ? "X" : name.substr(0, 1);
    }
    const std::u32string chars = core::decodeUtf8(sym);
    std::string first = core::encodeUtf8(chars.substr(0, 1));
    first = asciiUpper(first);
    return first + core::pythonLower(core::encodeUtf8(chars.substr(1)));
}

std::vector<std::string> columnOr(const CifLoop& loop, std::string_view item, std::size_t n) {
    std::vector<std::string> col = loop.column(item);
    if (col.empty()) {
        col.assign(n, "");
    }
    return col;
}

std::size_t heavyCount(const CifLoop& atoms) {
    const std::vector<std::string> names = atoms.column("atom_id");
    const std::vector<std::string> symbols = columnOr(atoms, "type_symbol", names.size());
    std::size_t n = 0;
    for (std::size_t k = 0; k < std::min(names.size(), symbols.size()); ++k) {
        const std::string el = elementOf(symbols[k], names[k]);
        if (el != "H" && el != "D") {
            ++n;
        }
    }
    return n;
}

std::size_t componentBlock(const CifDocument& doc, const std::string& sourceId, const Component* comp) {
    std::vector<std::size_t> blocks;
    for (std::size_t k = 0; k < doc.blocks.size(); ++k) {
        const CifLoop* atoms = doc.blocks[k].get("_chem_comp_atom");
        if (atoms && !atoms->rows.empty()) {
            blocks.push_back(k);
        }
    }
    if (blocks.empty()) {
        fail("no _chem_comp_atom loop found - is this a restraint dictionary?");
    }
    const auto names = [&] {
        std::vector<std::string> out;
        for (const std::size_t k : blocks) {
            out.push_back(doc.blocks[k].name);
        }
        return join(out, ", ");
    };
    if (!sourceId.empty()) {
        for (const std::size_t k : blocks) {
            const CifBlock& b = doc.blocks[k];
            if (core::pythonLower(b.name) == core::pythonLower("comp_" + sourceId) ||
                contains(b.get("_chem_comp_atom")->column("comp_id"), sourceId)) {
                return k;
            }
        }
        fail("no component " + repr(sourceId) + " in the dictionary (found: " + names() + ")");
    }
    if (blocks.size() == 1) {
        return blocks[0];
    }
    if (comp != nullptr) {  // the block whose heavy-atom count fits the reference
        const std::size_t full = comp->atoms.size(), residue = comp->atoms.size() - comp->leaving.size();
        std::vector<std::size_t> hits;
        for (const std::size_t k : blocks) {
            const std::size_t n = heavyCount(*doc.blocks[k].get("_chem_comp_atom"));
            if (n == full || n == residue) {
                hits.push_back(k);
            }
        }
        if (hits.size() == 1) {
            return hits[0];
        }
    }
    fail("the dictionary holds several components (" + names() + "); choose one with --source-id");
}

std::optional<int> formalCharge(const std::optional<std::string>& text) {
    if (!text) {
        return std::nullopt;
    }
    const std::optional<double> q = core::tryPythonFloat(*text);
    // DIC-01: a non-finite or out-of-range value is no formal charge (MDATK: OverflowError)
    if (!q || !std::isfinite(*q) || std::fabs(*q) > 1.0e6) {
        return std::nullopt;
    }
    const double rounded = std::nearbyint(*q);  // round(): ties to even
    if (std::fabs(*q - rounded) > 1e-6) {
        return std::nullopt;  // a partial charge, not a formal one
    }
    return static_cast<int>(rounded);
}

struct DictionaryGraph {
    MolGraph graph;
    Coords coords;
};

DictionaryGraph dictionaryGraph(const CifBlock& block) {
    const CifLoop& atoms = *block.get("_chem_comp_atom");
    const std::vector<std::string> names = atoms.column("atom_id");
    if (names.empty()) {
        fail("_chem_comp_atom has no atom_id column");
    }
    const std::vector<std::string> symbols = columnOr(atoms, "type_symbol", names.size());
    const std::vector<std::string> chargeColumn = atoms.column("charge");
    std::vector<std::optional<std::string>> chargesRaw(names.size());
    if (!chargeColumn.empty()) {
        chargesRaw.assign(chargeColumn.begin(), chargeColumn.end());
    }
    const std::size_t n = std::min(names.size(), symbols.size());
    std::vector<std::string> elements;
    for (std::size_t k = 0; k < n; ++k) {
        elements.push_back(elementOf(symbols[k], names[k]));
    }
    DictionaryGraph out;
    static constexpr std::array<std::array<std::string_view, 3>, 3> kCoords = {{
        {"x", "y", "z"},
        {"model_cartn_x", "model_cartn_y", "model_cartn_z"},
        {"pdbx_model_cartn_x_ideal", "pdbx_model_cartn_y_ideal", "pdbx_model_cartn_z_ideal"},
    }};
    for (const auto& cols : kCoords) {
        const std::array<std::vector<std::string>, 3> xyz = {atoms.column(cols[0]), atoms.column(cols[1]),
                                                             atoms.column(cols[2])};
        if (xyz[0].empty() || xyz[1].empty() || xyz[2].empty()) {
            continue;
        }
        Coords coords;
        const std::size_t m = std::min({names.size(), xyz[0].size(), xyz[1].size(), xyz[2].size()});
        for (std::size_t k = 0; k < m; ++k) {
            const auto x = core::tryPythonFloat(xyz[0][k]), y = core::tryPythonFloat(xyz[1][k]),
                       z = core::tryPythonFloat(xyz[2][k]);
            if (x && y && z) {
                coords.set(names[k], Vec{*x, *y, *z});
            }
        }
        const bool nonZero = std::any_of(coords.items().begin(), coords.items().end(), [](const auto& item) {
            return std::any_of(item.second.begin(), item.second.end(), [](double v) { return std::fabs(v) > 1e-6; });
        });
        if (!coords.empty() && nonZero) {
            out.coords = std::move(coords);
            break;
        }
    }
    {
        std::set<std::string> unique(names.begin(), names.end());
        if (unique.size() != names.size()) {
            fail("duplicate atom names in _chem_comp_atom");
        }
    }
    MolGraph& g = out.graph;
    std::map<std::string, int> index;
    OrderedMap<std::string> hydrogens;  // hydrogen -> its heavy atom ("" until bonded)
    const std::size_t rows = std::min(n, chargesRaw.size());
    for (std::size_t k = 0; k < rows; ++k) {
        const std::string& el = elements[k];
        if (el == "H" || el == "D") {
            hydrogens.set(names[k], "");
            continue;
        }
        index[names[k]] = g.addAtom(el, formalCharge(chargesRaw[k]), 0, names[k]);
    }
    const CifLoop* bonds = block.get("_chem_comp_bond");
    if (bonds == nullptr || bonds->rows.empty()) {
        fail("no _chem_comp_bond loop - cannot match atoms without bonds");
    }
    const std::vector<std::string> a1 = bonds->column("atom_id_1"), a2 = bonds->column("atom_id_2");
    std::vector<std::string> kinds = bonds->column("type");
    if (kinds.empty()) {
        kinds = bonds->column("value_order");
    }
    if (kinds.empty()) {
        kinds.assign(a1.size(), "");
    }
    const auto known = [&](const std::string& name) { return index.count(name) != 0 || hydrogens.has(name); };
    for (std::size_t k = 0; k < std::min({a1.size(), a2.size(), kinds.size()}); ++k) {
        const std::string& x = a1[k];
        const std::string& y = a2[k];
        if (!known(x) || !known(y)) {
            fail("bond " + x + "-" + y + " names an atom missing from _chem_comp_atom");
        }
        if (hydrogens.has(x) || hydrogens.has(y)) {
            const bool xh = hydrogens.has(x);
            const std::string& h = xh ? x : y;
            const std::string& heavy = xh ? y : x;
            if (hydrogens.has(heavy)) {
                fail("H-H bond " + x + "-" + y);
            }
            if (!hydrogens.get(h)->empty()) {
                fail("hydrogen " + h + " has two bonds");
            }
            hydrogens.set(h, heavy);
            const auto hi = static_cast<std::size_t>(index.at(heavy));
            ++g.hcount[hi];
            g.hnames[hi].push_back(h);
            continue;
        }
        static const std::map<std::string, int> kOrders = {
            {"single", 1}, {"sing", 1}, {"1", 1}, {"double", 2}, {"doub", 2}, {"2", 2},
            {"triple", 3}, {"trip", 3}, {"3", 3}};
        const auto it = kOrders.find(core::pythonLower(core::pythonStrip(kinds[k])));
        g.addBond(index.at(x), index.at(y), it == kOrders.end() ? 0 : it->second);
    }
    std::vector<std::string> lone;
    for (const auto& [h, heavy] : hydrogens.items()) {
        if (heavy.empty()) {
            lone.push_back(h);
        }
    }
    if (!lone.empty()) {
        fail("hydrogens without a bond: " + join(lone, ", "));
    }
    return out;
}

int hposCost(const std::string& kind, const Vec& h, const std::array<Vec, 3>& pts, bool mirrored) {
    if (kind == "+" || kind == "-") {
        return static_cast<int>((volume(pts[0], pts[1], pts[2], h) > 0) != ((kind == "+") != mirrored));
    }
    const bool cis = std::fabs(torsion(h, pts[0], pts[1], pts[2])) < 90.0;
    return static_cast<int>(cis != (kind == "cis"));
}

// _name_hydrogens: old hydrogen name -> reference name for every hydrogen of the dictionary
// (`mirrored`: D coordinates, so CH2 hydrogens are told apart by the mirrored rule).
void nameHydrogens(const Component& comp, const MolGraph& ref, const MolGraph& tgt, const std::vector<int>& m,
                   const Coords& coords, bool allowProtonation, std::vector<std::string>& notes, bool mirrored,
                   NameMap& out) {
    std::map<std::string, const HydrogenRule*> hpos;
    for (const HydrogenRule& r : comp.hpos) {
        hpos[r.hydrogen] = &r;  // a later rule of the same H wins
    }
    std::map<std::string, std::size_t> refIndex;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        refIndex[ref.names[i]] = i;  // the last index of a name wins
    }
    std::set<std::string> taken;
    for (const std::string& a : comp.residueAtoms()) {
        taken.insert(a);
    }
    for (const auto& a : comp.atoms) {
        taken.insert(a.name);
    }
    for (std::size_t i = 0; i < ref.size(); ++i) {
        const std::string& xname = ref.names[i];
        const std::vector<std::string>& gen = tgt.hnames[static_cast<std::size_t>(m[i])];
        const std::vector<std::string>& want = ref.hnames[i];
        if (gen.size() != want.size()) {
            if (!protonationName(xname) && !allowProtonation) {
                fail("the dictionary has " + std::to_string(gen.size()) + " H on " + xname + ", " + comp.code +
                     " has " + std::to_string(want.size()) +
                     " (different protonation / tautomer; --allow-protonation accepts it)");
            }
            notes.push_back(xname + " carries " + std::to_string(gen.size()) + " H (reference " +
                            std::to_string(want.size()) + ")");
        }
        std::vector<std::string> names(want.begin(), want.begin() + static_cast<std::ptrdiff_t>(
                                                                         std::min(want.size(), gen.size())));
        int k = 1;
        while (names.size() < gen.size()) {  // extra hydrogens: CCD-style name when there is one
            std::string base;
            if (const auto site = protonationName(xname)) {
                base = *site;
            } else {
                base = "H";
                for (const char c : xname) {
                    if (c != '\'' && c != '"') {
                        base += c;
                    }
                }
            }
            const std::string cand = k == 1 ? base : base + std::to_string(k);
            ++k;
            if (taken.count(cand) == 0) {
                names.push_back(cand);
                taken.insert(cand);
            }
        }
        const std::set<std::string> genSet(gen.begin(), gen.end()), wantSet(want.begin(), want.end());
        if (genSet == wantSet) {  // already the reference names (e.g. a CCD dictionary): keep them
            names = gen;
        } else if (gen.size() >= 2 && gen.size() == want.size() &&
                   std::all_of(want.begin(), want.end(), [&](const std::string& h) { return hpos.count(h) != 0; })) {
            const HydrogenRule& first = *hpos.at(want[0]);
            std::array<Vec, 3> pts{};
            bool havePts = true;
            const std::string* keys[] = {&first.x, &first.a, &first.b};
            for (std::size_t p = 0; p < 3; ++p) {
                const auto ri = refIndex.find(*keys[p]);
                const Vec* v = ri == refIndex.end()
                                   ? nullptr
                                   : coords.get(tgt.names[static_cast<std::size_t>(m[ri->second])]);
                if (v == nullptr) {
                    havePts = false;
                    break;
                }
                pts[p] = *v;
            }
            if (havePts && std::all_of(gen.begin(), gen.end(), [&](const std::string& h) { return coords.has(h); })) {
                std::vector<std::size_t> perm(want.size());
                std::iota(perm.begin(), perm.end(), std::size_t{0});
                std::optional<std::pair<int, std::vector<std::size_t>>> best;
                do {  // _permutations: itertools order
                    int cost = 0;
                    for (std::size_t j = 0; j < gen.size(); ++j) {
                        cost += hposCost(hpos.at(want[perm[j]])->kind, *coords.get(gen[j]), pts, mirrored);
                    }
                    if (!best || cost < best->first) {
                        best = std::make_pair(cost, perm);
                    }
                } while (std::next_permutation(perm.begin(), perm.end()));
                names.clear();
                for (const std::size_t j : best->second) {
                    names.push_back(want[j]);
                }
            }
        }
        for (std::size_t j = 0; j < std::min(gen.size(), names.size()); ++j) {
            out.set(gen[j], names[j]);
        }
    }
}

// _coordinate_hands: sugar centre -> +1 (L) / -1 (D) from the coordinates, in the centres' order.
std::vector<std::pair<std::string, int>> coordinateHands(const Component& comp, const Coords& coords,
                                                         const NameMap& mapping) {
    Coords renamed;
    for (const auto& [name, p] : coords.items()) {
        renamed.set(mapped(mapping, name), p);
    }
    std::vector<std::pair<std::string, int>> out;
    for (const ChiralityReference& c : comp.chirality()) {
        const Vec* pts[] = {renamed.get(c.centre), renamed.get(c.order[0]), renamed.get(c.order[1]),
                            renamed.get(c.order[2])};
        if (std::all_of(std::begin(pts), std::end(pts), [](const Vec* p) { return p != nullptr; })) {
            const double vol = volume(*pts[0], *pts[1], *pts[2], *pts[3]);
            if (std::fabs(vol) > 0.1) {
                out.emplace_back(c.centre, (vol > 0) == (c.sign > 0) ? 1 : -1);
            }
        }
    }
    return out;
}

void checkChirality(const Component& comp, const CifBlock& block, const std::vector<std::pair<std::string, int>>& hands,
                    const NameMap& mapping, std::vector<std::string>& notes) {
    const std::vector<ChiralityReference> centres = comp.chirality();
    const auto refOf = [&](const std::string& centre) -> const ChiralityReference* {
        for (const auto& c : centres) {
            if (c.centre == centre) {
                return &c;
            }
        }
        return nullptr;
    };
    std::set<std::string> checked;
    const CifLoop* chir = block.get("_chem_comp_chir");
    if (chir != nullptr && !chir->rows.empty()) {
        std::array<std::ptrdiff_t, 5> cols{};
        const std::string_view items[] = {"atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "volume_sign"};
        for (std::size_t k = 0; k < 5; ++k) {
            cols[k] = chir->col(items[k]);
        }
        if (*std::min_element(cols.begin(), cols.end()) >= 0) {
            for (const auto& row : chir->rows) {
                const auto cell = [&](std::size_t k) { return row.at(static_cast<std::size_t>(cols[k])); };
                const std::string centre = mapped(mapping, cell(0));
                const std::string b1 = mapped(mapping, cell(1)), b2 = mapped(mapping, cell(2)),
                                  b3 = mapped(mapping, cell(3));
                const std::string vs = core::pythonLower(cell(4));
                const ChiralityReference* ref = refOf(centre);
                if (ref == nullptr || startsWith(vs, "both")) {
                    continue;
                }
                const std::vector<std::string> order = orderNames(*ref);
                std::vector<std::string> rest;
                for (const std::string& name : order) {
                    if (name != b1 && name != b2 && name != b3) {
                        rest.push_back(name);
                    }
                }
                const int rowSign = startsWith(vs, "pos") ? 1 : startsWith(vs, "neg") ? -1 : 0;
                if (rest.size() != 1 || rowSign == 0) {
                    continue;
                }
                if (expectedSign(order, ref->sign, {b1, b2, b3, rest[0]}) != rowSign) {
                    fail("_chem_comp_chir: " + centre + " has the opposite configuration (" + cell(4) + "); " +
                         comp.code + " needs an L sugar - check the SMILES");
                }
                checked.insert(centre);
            }
        }
    }
    std::vector<std::string> wrong;
    for (const auto& [centre, hand] : hands) {
        if (hand < 0) {
            wrong.push_back(centre);
        }
    }
    if (!wrong.empty() && checked.empty()) {
        fail("the coordinates have the opposite configuration at " + join(wrong, ", ") + "; " + comp.code +
             " needs an L sugar - check the SMILES");
    }
    if (!wrong.empty()) {
        notes.push_back("WARNING: the ideal coordinates have the D configuration at " + join(wrong, ", ") +
                        " although the chirality restraints are L (programs that build the residue from these "
                        "coordinates place a D sugar)");
    }
    std::vector<std::string> missing;
    for (const auto& c : centres) {
        const bool inHands = std::any_of(hands.begin(), hands.end(), [&](const auto& h) { return h.first == c.centre; });
        if (checked.count(c.centre) == 0 && !inHands) {
            missing.push_back(c.centre);
        }
    }
    if (!missing.empty()) {
        notes.push_back("configuration not checked (no chirality rows or coordinates) at " + join(missing, ", "));
    }
}

std::string programOf(const CifDocument& doc) {
    std::string text = core::pythonLower(join(doc.header, " "));
    for (const CifBlock& b : doc.blocks) {
        if (const CifLoop* gen = b.get("_pdbx_chem_comp_description_generator")) {
            text += " " + core::pythonLower(join(gen->column("program_name"), " "));
        }
    }
    if (text.find("elbow") != std::string::npos) {
        return "phenix.elbow";
    }
    if (text.find("acedrg") != std::string::npos) {
        return "acedrg";
    }
    return "";
}

std::string_view itemName(std::string_view tag) {
    const std::size_t dot = tag.find('.');
    return dot == std::string_view::npos ? tag : tag.substr(dot + 1);
}

// _NAME_ITEMS: atom_id | alt_atom_id | pdbx_component_atom_id | atom_id_<digits> | atom_id_centre
// (case-insensitive)
bool isNameItem(std::string_view tag) {
    const std::string item = core::pythonLower(itemName(tag));
    if (item == "atom_id" || item == "alt_atom_id" || item == "pdbx_component_atom_id" || item == "atom_id_centre") {
        return true;
    }
    if (!startsWith(item, "atom_id_") || item.size() == 8) {
        return false;
    }
    const std::u32string digits = core::decodeUtf8(std::string_view(item).substr(8));
    return std::all_of(digits.begin(), digits.end(), core::isPythonDigit);
}

// _COMP_ITEMS: comp_id | pdbx_component_comp_id (case-insensitive)
bool isCompItem(std::string_view tag) {
    const std::string item = core::pythonLower(itemName(tag));
    return item == "comp_id" || item == "pdbx_component_comp_id";
}

// _OUR_HEADER: a header line MDATK (or MDBridge) wrote
bool isOurHeader(const std::string& line) {
    return line.find("renamed by MDATK") != std::string::npos || startsWith(line, "# _chem_comp.group ");
}

bool isAsciiAlnum(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
}

}  // namespace

FinishResult finishDictionary(std::string_view text, std::string_view codeIn, const FinishOptions& options) {
    // PDB residue names are upper case; Coot / Phenix match ids exactly
    const std::string code = asciiUpper(core::pythonStrip(codeIn));
    if (code.empty() || code.size() > 5 || !std::all_of(code.begin(), code.end(), isAsciiAlnum)) {
        fail("residue code " + repr(code) + ": use 1-5 letters / digits");
    }
    const Component* comp = options.reference != nullptr ? options.reference : findComponent(code);
    if (options.keepNames) {
        comp = nullptr;
    }
    CifDocument doc = io::parseCif(text);
    const std::size_t blockIndex = componentBlock(doc, options.sourceId, comp);
    CifBlock* block = &doc.blocks[blockIndex];
    std::string oldId;
    {
        const std::vector<std::string> ids = block->get("_chem_comp_atom")->column("comp_id");
        oldId = ids.empty() ? std::string() : ids[0];
        if (oldId.empty()) {
            oldId = block->name;
            if (oldId.size() >= 5 && core::pythonLower(oldId.substr(0, 5)) == "comp_") {
                oldId = oldId.substr(5);
            }
        }
    }
    std::vector<std::string> notes;
    const std::string program = programOf(doc);
    NameMap mapping;
    std::set<std::string> removed;
    if (comp != nullptr) {
        DictionaryGraph dg = dictionaryGraph(*block);
        const MolGraph& tgt = dg.graph;
        const bool withLeaving = tgt.size() != comp->atoms.size() - comp->leaving.size() || comp->leaving.empty();
        const MolGraph ref = comp->graph(withLeaving);
        const std::optional<std::vector<int>> m = matchGraph(ref, tgt);
        if (!m) {
            fail(describeMismatch(*comp, tgt, "the dictionary " + oldId));
        }
        for (std::size_t i = 0; i < ref.size(); ++i) {
            mapping.set(tgt.names[static_cast<std::size_t>((*m)[i])], ref.names[i]);
        }
        const auto hands = coordinateHands(*comp, dg.coords, mapping);
        const bool mirrored = !hands.empty() && std::all_of(hands.begin(), hands.end(),
                                                            [](const auto& h) { return h.second < 0; });
        nameHydrogens(*comp, ref, tgt, *m, dg.coords, options.allowProtonation, notes, mirrored, mapping);
        checkChirality(*comp, *block, hands, mapping, notes);
        for (const auto& [o, n] : mapping.items()) {
            if (comp->isLeaving(n)) {
                removed.insert(o);
            }
        }
        for (std::size_t i = 0; i < ref.size(); ++i) {
            if (comp->isLeaving(ref.names[i])) {
                const auto& hs = tgt.hnames[static_cast<std::size_t>((*m)[i])];
                removed.insert(hs.begin(), hs.end());
            }
        }
        const auto changed = std::count_if(mapping.items().begin(), mapping.items().end(),
                                           [](const auto& item) { return item.first != item.second; });
        std::string note = std::to_string(changed) + " of " + std::to_string(mapping.size()) +
                           " atom names changed to the " + comp->parent + " (CCD) names";
        if (!removed.empty()) {
            std::vector<std::string> gone;
            for (const std::string& o : removed) {
                gone.push_back(*mapping.get(o));
            }
            std::sort(gone.begin(), gone.end());
            note += "; removed leaving atoms " + join(gone, ", ");
        }
        notes.push_back(note);
        std::set<std::string> values;
        for (const auto& item : mapping.items()) {
            values.insert(item.second);
        }
        if (values.size() != mapping.size()) {
            fail("internal error: two atoms got the same name");
        }
    } else if (!options.keepNames) {
        notes.push_back(code + " has no reference atom names; the program's names are kept");
    }
    NameMap newNames;
    for (const auto& [o, n] : mapping.items()) {
        if (removed.count(o) == 0) {
            newNames.set(o, n);
        }
    }
    // rename / delete inside the component block
    for (CifLoop& item : block->items) {
        std::vector<std::size_t> nameCols, compCols;
        for (std::size_t k = 0; k < item.tags.size(); ++k) {
            if (isNameItem(item.tags[k])) {
                nameCols.push_back(k);
            }
            if (isCompItem(item.tags[k])) {
                compCols.push_back(k);
            }
        }
        std::vector<std::vector<std::string>> rows;
        for (std::vector<std::string>& row : item.rows) {
            if (std::any_of(nameCols.begin(), nameCols.end(), [&](std::size_t k) { return removed.count(row[k]) != 0; })) {
                continue;
            }
            for (const std::size_t k : nameCols) {
                row[k] = mapped(newNames, row[k]);
            }
            for (const std::size_t k : compCols) {
                if (row[k] == oldId || row[k].empty() || row[k] == "." || row[k] == "?" || row[k] == ".?") {
                    row[k] = code;
                }
            }
            rows.push_back(std::move(row));
        }
        item.rows = std::move(rows);
    }
    std::erase_if(block->items, [](const CifLoop& it) { return it.rows.empty() && it.pairs; });
    block->name = "comp_" + code;
    const CifLoop* atomsLoop = block->get("_chem_comp_atom");
    if (atomsLoop == nullptr) {
        fail("internal error: the component lost its _chem_comp_atom loop");
    }
    const std::size_t nAll = atomsLoop->rows.size();
    const std::size_t nHeavy = heavyCount(*atomsLoop);
    const std::vector<std::string> got = atomsLoop->column("atom_id");
    std::string title = options.name;
    if (title.empty()) {
        if (comp != nullptr) {
            title = comp->name;
        } else if (const Component* known = findComponent(code)) {
            title = known->name;
        } else {
            title = code;
        }
    }
    // the component list
    CifLoop* clist = nullptr;
    for (CifBlock& b : doc.blocks) {
        CifLoop* loop = b.get("_chem_comp");
        if (loop != nullptr && loop->col("id") >= 0) {
            clist = loop;
            break;
        }
    }
    if (clist == nullptr) {
        CifBlock list;
        list.name = "comp_list";
        CifLoop loop;
        loop.tags = {"_chem_comp.id", "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
                     "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh", "_chem_comp.desc_level"};
        loop.rows = {{oldId, oldId, "", "", "", "", "."}};
        list.items.push_back(std::move(loop));
        doc.blocks.insert(doc.blocks.begin(), std::move(list));
        block = nullptr;  // the vector moved
        clist = &doc.blocks.front().items.front();
    }
    const auto idCol = static_cast<std::size_t>(clist->col("id"));
    std::size_t rowIndex = clist->rows.size();
    for (std::size_t r = 0; r < clist->rows.size(); ++r) {
        const std::string& id = clist->rows[r].at(idCol);
        if (id == oldId || id == code) {
            rowIndex = r;
            break;
        }
    }
    if (rowIndex == clist->rows.size()) {
        clist->rows.emplace_back(clist->tags.size(), "");
    }
    const std::pair<std::string_view, std::string> values[] = {
        {"id", code},          {"three_letter_code", code.substr(0, 3)},       {"name", title},
        {"group", options.group}, {"number_atoms_all", std::to_string(nAll)}, {"number_atoms_nh", std::to_string(nHeavy)}};
    for (const auto& [item, value] : values) {
        std::ptrdiff_t k = clist->col(item);
        if (k < 0) {
            clist->tags.push_back("_chem_comp." + std::string(item));
            for (auto& r : clist->rows) {
                r.emplace_back(".");
            }
            k = static_cast<std::ptrdiff_t>(clist->tags.size()) - 1;
        }
        clist->rows[rowIndex].at(static_cast<std::size_t>(k)) = value;
    }
    for (auto& r : clist->rows) {
        for (std::string& v : r) {
            if (v.empty()) {
                v = ".";
            }
        }
    }
    const std::vector<std::string> oldHeader = doc.header;
    doc.header = {"# " + code + ": restraint dictionary" + (program.empty() ? "" : " from " + program) +
                      ", renamed by MDATK dictionaries.py" +
                      (comp != nullptr ? " to the wwPDB CCD atom names of " + comp->parent : ""),
                  "# _chem_comp.group " + options.group};
    for (const std::string& h : oldHeader) {
        if (core::pythonStrip(h) != "#" && !isOurHeader(h)) {
            doc.header.push_back(h);
        }
    }
    FinishResult result;
    result.text = io::writeCif(doc);
    if (comp != nullptr) {
        std::vector<std::string> missing;
        for (const std::string& a : comp->residueAtoms()) {
            if (!contains(got, a)) {
                missing.push_back(a);
            }
        }
        if (!missing.empty()) {
            notes.push_back("atoms of the model without restraints: " + join(missing, ", "));
        }
    }
    for (const auto& [o, n] : mapping.items()) {
        result.mapping.emplace_back(o, removed.count(o) != 0 ? std::string() : n);
    }
    result.notes = std::move(notes);
    result.sourceId = oldId;
    result.program = program;
    return result;
}

std::string renameDictionary(std::string_view text, std::string_view code, const FinishOptions& options) {
    return finishDictionary(text, code, options).text;
}

}  // namespace mdbridge::chem
