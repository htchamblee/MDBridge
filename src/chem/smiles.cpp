// Molecular graphs, SMILES and graph matching of MDATK's mdatk/dictionaries.py (MDATK's own code,
// ported): parse_smiles, MolGraph, _all_matches, _score, match_graph.

#include "mdbridge/chem/smiles.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <stdexcept>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::chem {

using core::decodeUtf8;
using core::encodeUtf8;
using core::isPythonDigit;
using core::isPythonSpace;

// --- MolGraph ------------------------------------------------------------------------------------

int MolGraph::addAtom(std::string element, std::optional<int> charge, int hydrogens, std::string name,
                      std::vector<std::string> hydrogenNames) {
    elements.push_back(std::move(element));
    charges.push_back(charge);
    hcount.push_back(hydrogens);
    names.push_back(std::move(name));
    hnames.push_back(std::move(hydrogenNames));
    adj.emplace_back();
    return static_cast<int>(elements.size()) - 1;
}

namespace {

void setNeighbour(std::vector<std::pair<int, int>>& list, int other, int order) {
    for (auto& entry : list) {
        if (entry.first == other) {
            entry.second = order;
            return;
        }
    }
    list.emplace_back(other, order);
}

bool hasNeighbour(const std::vector<std::pair<int, int>>& list, int other) {
    return std::any_of(list.begin(), list.end(), [&](const auto& entry) { return entry.first == other; });
}

}  // namespace

void MolGraph::addBond(int a, int b, int order) {
    setNeighbour(adj[static_cast<std::size_t>(a)], b, order);
    setNeighbour(adj[static_cast<std::size_t>(b)], a, order);
}

void MolGraph::removeBond(int a, int b) {
    auto drop = [](std::vector<std::pair<int, int>>& list, int other) {
        list.erase(std::remove_if(list.begin(), list.end(), [&](const auto& e) { return e.first == other; }),
                   list.end());
    };
    drop(adj[static_cast<std::size_t>(a)], b);
    drop(adj[static_cast<std::size_t>(b)], a);
}

std::optional<int> MolGraph::bond(int a, int b) const {
    for (const auto& [other, order] : adj[static_cast<std::size_t>(a)]) {
        if (other == b) {
            return order;
        }
    }
    return std::nullopt;
}

std::string MolGraph::formula() const {
    std::map<std::string, long long> counts;  // sorted by code point, as Python's sorted(str)
    for (const std::string& el : elements) {
        ++counts[el];
    }
    long long h = 0;
    for (const int n : hcount) {
        h += n;
    }
    if (h != 0) {
        counts["H"] += h;
    }
    std::vector<std::string> keys;
    if (counts.count("C") != 0) {
        keys.emplace_back("C");
        if (counts.count("H") != 0) {
            keys.emplace_back("H");
        }
    }
    for (const auto& [key, _n] : counts) {
        if (std::find(keys.begin(), keys.end(), key) == keys.end()) {
            keys.push_back(key);
        }
    }
    std::string text;
    for (const std::string& key : keys) {
        text += key;
        if (counts[key] > 1) {
            text += std::to_string(counts[key]);
        }
    }
    const bool known = std::all_of(charges.begin(), charges.end(), [](const auto& q) { return q.has_value(); });
    if (known) {
        long long q = 0;
        for (const auto& c : charges) {
            q += *c;
        }
        if (q != 0) {
            text += q > 0 ? "+" : "-";
            if (std::llabs(q) > 1) {
                text += std::to_string(std::llabs(q));
            }
        }
    }
    return text;
}

// --- SMILES ----------------------------------------------------------------------------------------

namespace {

bool isUpper(char32_t c) { return c >= U'A' && c <= U'Z'; }
bool isLower(char32_t c) { return c >= U'a' && c <= U'z'; }

std::string utf8(std::u32string_view text) { return encodeUtf8(text); }

[[noreturn]] void fail(const std::string& message) { throw DictionaryError(message); }

// The length of the _SMILES_TOKEN match at `pos` (0 if none):
// \[[^\]]*\]|Br|Cl|[BCNOPSFI]|[bcnops]|%\d\d|\d|[-=#$:/\\]|[().]|\*
std::size_t tokenAt(std::u32string_view text, std::size_t pos) {
    const char32_t c = text[pos];
    const auto next = [&](std::size_t k) -> char32_t { return pos + k < text.size() ? text[pos + k] : U'\0'; };
    if (c == U'[') {
        const std::size_t close = text.find(U']', pos + 1);
        return close == std::u32string_view::npos ? 0 : close - pos + 1;
    }
    if ((c == U'B' && next(1) == U'r') || (c == U'C' && next(1) == U'l')) {
        return 2;
    }
    if (std::u32string_view(U"BCNOPSFIbcnops-=#$:/\\().*").find(c) != std::u32string_view::npos) {
        return 1;
    }
    if (c == U'%') {
        return (pos + 2 < text.size() && isPythonDigit(text[pos + 1]) && isPythonDigit(text[pos + 2])) ? 3 : 0;
    }
    return isPythonDigit(c) ? 1 : 0;
}

int bondOrder(char32_t symbol) {
    switch (symbol) {
    case U'=': return 2;
    case U'#': return 3;
    case U'$': return 4;
    case U':': return 0;
    default: return 1;  // - / backslash
    }
}

bool isBondSymbol(std::u32string_view tok) {
    return tok.size() == 1 && std::u32string_view(U"-=#$:/\\").find(tok[0]) != std::u32string_view::npos;
}

// int() of a run of (Unicode) decimal digits; values beyond what a molecule can hold are refused
// as an unreadable atom (MDBridge; Python has unbounded integers).
long long digitValue(std::u32string_view digits, std::u32string_view token) {
    const std::int64_t v = [&] {
        try {
            return core::pythonInt(utf8(digits));
        } catch (const std::overflow_error&) {
            return std::int64_t{-1};
        }
    }();
    if (v < 0 || v > 1000000) {
        fail("SMILES: cannot read the atom " + utf8(token));
    }
    return v;
}

struct BracketAtom {
    std::u32string symbol;
    int hydrogens = 0;
    int charge = 0;
    std::string stereo;
};

// _BRACKET: ^\[(\d*)([A-Z][a-z]?|[a-z][a-z]?|\*)(@@|@(?:TH[12])?)?(H\d*)?([+-]+\d*)?(?::\d+)?\]$
// (each part greedy; no later part can start with what an earlier one leaves, so no backtracking
// changes the match).
std::optional<BracketAtom> bracketAtom(std::u32string_view tok) {
    const std::u32string_view body = tok.substr(1, tok.size() - 2);
    std::size_t i = 0;
    const auto at = [&](std::size_t k) -> char32_t { return k < body.size() ? body[k] : U'\0'; };
    while (i < body.size() && isPythonDigit(body[i])) {
        ++i;  // isotope
    }
    BracketAtom atom;
    if (isUpper(at(i)) || isLower(at(i))) {  // [A-Z][a-z]? | [a-z][a-z]?
        atom.symbol = body.substr(i, 1);
        ++i;
        if (isLower(at(i))) {
            atom.symbol += body[i];
            ++i;
        }
    } else if (at(i) == U'*') {
        atom.symbol = U"*";
        ++i;
    } else {
        return std::nullopt;
    }
    if (at(i) == U'@') {
        if (at(i + 1) == U'@') {
            atom.stereo = "@@";
            i += 2;
        } else if (at(i + 1) == U'T' && at(i + 2) == U'H' && (at(i + 3) == U'1' || at(i + 3) == U'2')) {
            atom.stereo = at(i + 3) == U'2' ? "@@" : "@";
            i += 4;
        } else {
            atom.stereo = "@";
            i += 1;
        }
    }
    if (at(i) == U'H') {
        ++i;
        const std::size_t start = i;
        while (i < body.size() && isPythonDigit(body[i])) {
            ++i;
        }
        atom.hydrogens = i > start ? static_cast<int>(digitValue(body.substr(start, i - start), tok)) : 1;
    }
    if (at(i) == U'+' || at(i) == U'-') {
        const std::size_t start = i;
        while (at(i) == U'+' || at(i) == U'-') {
            ++i;
        }
        const std::size_t signs = i - start;
        const std::size_t digitsStart = i;
        while (i < body.size() && isPythonDigit(body[i])) {
            ++i;
        }
        const long long magnitude = i > digitsStart
                                        ? digitValue(body.substr(digitsStart, i - digitsStart), tok)
                                        : static_cast<long long>(signs);
        atom.charge = static_cast<int>(body[start] == U'+' ? magnitude : -magnitude);
    }
    if (at(i) == U':') {
        const std::size_t start = ++i;
        while (i < body.size() && isPythonDigit(body[i])) {
            ++i;
        }
        if (i == start) {
            return std::nullopt;
        }
    }
    if (i != body.size()) {
        return std::nullopt;
    }
    return atom;
}

// The valences of the organic subset (_VALENCE); empty for other elements.
std::vector<int> valences(const std::string& el) {
    static const std::map<std::string, std::vector<int>> table = {
        {"B", {3}},  {"C", {4}},  {"N", {3, 5}}, {"O", {2}},  {"P", {3, 5}},
        {"S", {2, 4, 6}}, {"F", {1}}, {"Cl", {1}}, {"Br", {1}}, {"I", {1}}};
    const auto it = table.find(el);
    return it == table.end() ? std::vector<int>{} : it->second;
}

// A neighbour of an atom while parsing: an atom index (>= 0) or the placeholder of an open ring
// bond (-2 - id of that opening).
constexpr int placeholder(int openId) { return -2 - openId; }

}  // namespace

MolGraph parseSmiles(std::string_view smiles) {
    const std::u32string all = decodeUtf8(smiles);
    // (smiles or "").strip().split()[0]
    std::size_t b = 0;
    while (b < all.size() && isPythonSpace(all[b])) {
        ++b;
    }
    std::size_t e = b;
    while (e < all.size() && !isPythonSpace(all[e])) {
        ++e;
    }
    const std::u32string text = all.substr(b, e - b);
    if (text.empty()) {
        fail("empty SMILES");
    }
    std::vector<std::u32string> tokens;
    std::size_t pos = 0;
    while (pos < text.size()) {
        const std::size_t len = tokenAt(text, pos);
        if (len == 0) {
            std::size_t next = pos + 1;
            while (next < text.size() && tokenAt(text, next) == 0) {
                ++next;
            }
            fail("SMILES: cannot read " + core::pythonStringRepr(utf8(text.substr(pos, next - pos))) +
                 " at position " + std::to_string(pos + 1));
        }
        tokens.push_back(text.substr(pos, len));
        pos += len;
    }

    std::vector<std::string> el;
    std::vector<bool> arom;
    std::vector<int> charge;
    std::vector<std::optional<int>> hexp;  // bracket H count, nullopt = implicit
    std::vector<std::string> tag;
    std::vector<std::vector<int>> nbrs;
    std::vector<bool> hasPrev;
    std::vector<std::pair<std::pair<int, int>, int>> bonds;  // insertion order
    std::map<std::pair<int, int>, std::size_t> bondIndex;
    int prev = -1;
    std::vector<int> stack;
    struct OpenRing {
        int atom;
        std::optional<char32_t> symbol;
        int id;
    };
    std::map<std::u32string, OpenRing> rings;
    int openings = 0;
    std::optional<char32_t> pending;

    auto bond = [&](int x, int y, std::optional<char32_t> symbol) {
        const std::pair<int, int> key{std::min(x, y), std::max(x, y)};
        if (bondIndex.count(key) != 0 || x == y) {
            fail("SMILES: duplicate bond / ring closure");
        }
        int order = 1;
        if (!symbol) {
            order = (arom[static_cast<std::size_t>(x)] && arom[static_cast<std::size_t>(y)]) ? 0 : 1;
        } else {
            order = bondOrder(*symbol);
        }
        bondIndex[key] = bonds.size();
        bonds.push_back({key, order});
    };

    for (const std::u32string& tok : tokens) {
        if (tok == U"(") {
            stack.push_back(prev);
        } else if (tok == U")") {
            if (stack.empty()) {
                fail("SMILES: unbalanced ')'");
            }
            prev = stack.back();
            stack.pop_back();
        } else if (tok == U".") {
            prev = -1;
        } else if (isBondSymbol(tok)) {
            pending = tok[0];
        } else if (tok[0] == U'%' || isPythonDigit(tok[0])) {
            if (prev < 0) {
                fail("SMILES: ring closure without an atom");
            }
            const std::u32string digit = tok[0] == U'%' ? tok.substr(1) : tok;
            const auto it = rings.find(digit);
            if (it != rings.end()) {
                const OpenRing open = it->second;
                rings.erase(it);
                bond(prev, open.atom, pending ? pending : open.symbol);
                nbrs[static_cast<std::size_t>(prev)].push_back(open.atom);
                auto& other = nbrs[static_cast<std::size_t>(open.atom)];
                *std::find(other.begin(), other.end(), placeholder(open.id)) = prev;
            } else {
                rings[digit] = OpenRing{prev, pending, openings};
                nbrs[static_cast<std::size_t>(prev)].push_back(placeholder(openings));
                ++openings;
            }
            pending.reset();
        } else {
            std::u32string sym;
            std::optional<int> h;
            int q = 0;
            std::string stereo;
            if (tok[0] == U'[') {
                const std::optional<BracketAtom> atom = bracketAtom(tok);
                if (!atom) {
                    fail("SMILES: cannot read the atom " + utf8(tok));
                }
                sym = atom->symbol;
                h = atom->hydrogens;
                q = atom->charge;
                stereo = atom->stereo;
            } else {
                sym = tok;
            }
            const bool aromatic = isLower(sym[0]);
            std::string element = utf8(sym);
            if (sym != U"*") {  // str.capitalize() (ASCII symbols)
                for (std::size_t k = 0; k < element.size(); ++k) {
                    const char c = element[k];
                    element[k] = k == 0 ? static_cast<char>(c >= 'a' && c <= 'z' ? c - 32 : c)
                                        : static_cast<char>(c >= 'A' && c <= 'Z' ? c + 32 : c);
                }
            }
            const int idx = static_cast<int>(el.size());
            el.push_back(element);
            arom.push_back(aromatic);
            charge.push_back(q);
            hexp.push_back(h);
            tag.push_back(stereo);
            nbrs.emplace_back();
            hasPrev.push_back(prev >= 0);
            if (prev >= 0) {
                bond(prev, idx, pending);
                nbrs[static_cast<std::size_t>(prev)].push_back(idx);
                nbrs[static_cast<std::size_t>(idx)].push_back(prev);
            }
            pending.reset();
            prev = idx;
        }
    }
    if (!rings.empty()) {
        std::string list;
        for (const auto& [digit, _open] : rings) {  // sorted by code point
            list += (list.empty() ? "" : ", ") + utf8(digit);
        }
        fail("SMILES: unclosed ring bond " + list);
    }
    if (!stack.empty()) {
        fail("SMILES: unbalanced '('");
    }

    // implicit hydrogens of organic-subset atoms
    const std::size_t n = el.size();
    std::vector<int> totalH(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (hexp[i]) {
            totalH[i] = *hexp[i];
            continue;
        }
        int used = 0;
        for (const int j : nbrs[i]) {
            const int order = bonds[bondIndex.at({std::min(static_cast<int>(i), j), std::max(static_cast<int>(i), j)})]
                                  .second;
            used += order != 0 ? order : 1;
        }
        const std::vector<int> vals = valences(el[i]);
        if (vals.empty()) {
            fail("SMILES: " + el[i] + " must be written in brackets");
        }
        if (arom[i]) {
            totalH[i] = std::max(0, vals[0] - used - 1);
        } else {
            const auto fit = std::find_if(vals.begin(), vals.end(), [&](int v) { return v >= used; });
            totalH[i] = fit != vals.end() ? *fit - used : 0;
        }
    }
    // fold explicit hydrogen atoms ([H]) into their neighbour
    std::vector<int> renumber(n, -1);
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < n; ++i) {
        if (!(el[i] == "H" && nbrs[i].size() == 1)) {
            renumber[i] = static_cast<int>(keep.size());
            keep.push_back(i);
        }
    }
    MolGraph g;
    for (const std::size_t i : keep) {
        g.addAtom(el[i], charge[i], totalH[i]);
    }
    for (const auto& [key, order] : bonds) {
        const int na = renumber[static_cast<std::size_t>(key.first)];
        const int nb = renumber[static_cast<std::size_t>(key.second)];
        if (na >= 0 && nb >= 0) {
            g.addBond(na, nb, order);
        } else {
            const int heavy = na >= 0 ? na : nb;
            if (heavy >= 0) {
                ++g.hcount[static_cast<std::size_t>(heavy)];
            }
        }
    }
    for (const std::size_t i : keep) {
        if (tag[i].empty()) {
            continue;
        }
        std::vector<int> order;
        for (const int j : nbrs[i]) {
            const int r = renumber[static_cast<std::size_t>(j)];
            order.push_back(r >= 0 ? r : MolGraph::kHydrogen);
        }
        if (hexp[i] && *hexp[i] > 0) {  // the bracket H follows the preceding atom (or comes first)
            for (int k = 0; k < *hexp[i]; ++k) {
                const std::size_t at = hasPrev[i] ? std::min<std::size_t>(1, order.size()) : 0;
                order.insert(order.begin() + static_cast<std::ptrdiff_t>(at), MolGraph::kHydrogen);
            }
        }
        if (order.size() == 4) {
            g.stereo[renumber[i]] = MolGraph::Stereo{order, tag[i]};
        }
    }
    return g;
}

// --- graph matching --------------------------------------------------------------------------------

namespace {

using Label = std::pair<std::string, std::size_t>;  // (element, degree)

class Matcher {
public:
    Matcher(const MolGraph& ref, const MolGraph& tgt) : m_ref(ref), m_tgt(tgt) {}

    // Calls visit(mapping) for each isomorphism in MDATK's order, at most kMaxMatches times.
    template <typename Visit>
    void run(Visit&& visit) {
        const std::size_t n = m_ref.size();
        if (n != m_tgt.size()) {
            return;
        }
        std::vector<std::string> er = m_ref.elements, et = m_tgt.elements;
        std::sort(er.begin(), er.end());
        std::sort(et.begin(), et.end());
        if (er != et) {
            return;
        }
        std::size_t degR = 0, degT = 0;
        for (std::size_t i = 0; i < n; ++i) {
            degR += m_ref.adj[i].size();
            degT += m_tgt.adj[i].size();
        }
        if (degR != degT) {
            return;
        }
        for (std::size_t i = 0; i < n; ++i) {
            m_labelR.emplace_back(m_ref.elements[i], m_ref.adj[i].size());
            m_labelT.emplace_back(m_tgt.elements[i], m_tgt.adj[i].size());
        }
        std::vector<Label> lr = m_labelR, lt = m_labelT;
        std::sort(lr.begin(), lr.end());
        std::sort(lt.begin(), lt.end());
        if (lr != lt) {
            return;
        }
        std::map<Label, std::size_t> freq;
        for (const Label& lab : m_labelT) {
            ++freq[lab];
        }
        // visiting order: rarest label first, then always a neighbour of the mapped part (BFS)
        const auto key = [&](int i) {
            return std::make_pair(freq[m_labelR[static_cast<std::size_t>(i)]], i);
        };
        std::vector<bool> seen(n, false);
        while (m_order.size() < n) {
            int start = -1;
            for (int i = 0; i < static_cast<int>(n); ++i) {
                if (!seen[static_cast<std::size_t>(i)] && (start < 0 || key(i) < key(start))) {
                    start = i;
                }
            }
            std::deque<int> queue{start};
            seen[static_cast<std::size_t>(start)] = true;
            while (!queue.empty()) {
                const int i = queue.front();
                queue.pop_front();
                m_order.push_back(i);
                std::vector<int> nb;
                for (const auto& [j, _o] : m_ref.adj[static_cast<std::size_t>(i)]) {
                    nb.push_back(j);
                }
                std::sort(nb.begin(), nb.end(), [&](int x, int y) { return key(x) < key(y); });
                for (const int j : nb) {
                    if (!seen[static_cast<std::size_t>(j)]) {
                        seen[static_cast<std::size_t>(j)] = true;
                        queue.push_back(j);
                    }
                }
            }
        }
        m_mapping.assign(n, -1);
        m_used.assign(n, false);
        extend(0, visit);
    }

private:
    template <typename Visit>
    void extend(std::size_t k, Visit& visit) {
        if (m_found >= kMaxMatches) {
            return;
        }
        const std::size_t n = m_ref.size();
        if (k == n) {
            ++m_found;
            visit(m_mapping);
            return;
        }
        const int r = m_order[k];
        std::vector<int> mappedNbrs;
        for (const auto& [j, _o] : m_ref.adj[static_cast<std::size_t>(r)]) {
            if (m_mapping[static_cast<std::size_t>(j)] >= 0) {
                mappedNbrs.push_back(j);
            }
        }
        std::vector<int> cands;
        if (!mappedNbrs.empty()) {
            const int anchor = m_mapping[static_cast<std::size_t>(mappedNbrs[0])];
            for (const auto& [t, _o] : m_tgt.adj[static_cast<std::size_t>(anchor)]) {
                if (!m_used[static_cast<std::size_t>(t)]) {
                    cands.push_back(t);
                }
            }
        } else {
            for (int t = 0; t < static_cast<int>(n); ++t) {
                if (!m_used[static_cast<std::size_t>(t)]) {
                    cands.push_back(t);
                }
            }
        }
        std::sort(cands.begin(), cands.end());
        for (const int t : cands) {
            const auto ts = static_cast<std::size_t>(t);
            if (m_labelT[ts] != m_labelR[static_cast<std::size_t>(r)]) {
                continue;
            }
            const auto& tadj = m_tgt.adj[ts];
            if (std::any_of(mappedNbrs.begin(), mappedNbrs.end(),
                            [&](int j) { return !hasNeighbour(tadj, m_mapping[static_cast<std::size_t>(j)]); })) {
                continue;
            }
            const auto usedNbrs = static_cast<std::size_t>(std::count_if(
                tadj.begin(), tadj.end(), [&](const auto& e) { return m_used[static_cast<std::size_t>(e.first)]; }));
            if (usedNbrs != mappedNbrs.size()) {
                continue;
            }
            m_mapping[static_cast<std::size_t>(r)] = t;
            m_used[ts] = true;
            extend(k + 1, visit);
            m_mapping[static_cast<std::size_t>(r)] = -1;
            m_used[ts] = false;
        }
    }

    const MolGraph& m_ref;
    const MolGraph& m_tgt;
    std::vector<Label> m_labelR, m_labelT;
    std::vector<int> m_order;
    std::vector<int> m_mapping;
    std::vector<bool> m_used;
    std::size_t m_found = 0;
};

}  // namespace

std::vector<std::vector<int>> allMatches(const MolGraph& ref, const MolGraph& tgt) {
    std::vector<std::vector<int>> out;
    Matcher(ref, tgt).run([&](const std::vector<int>& m) { out.push_back(m); });
    return out;
}

int matchScore(const MolGraph& ref, const MolGraph& tgt, const std::vector<int>& m) {
    int bad = 0;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        const auto t = static_cast<std::size_t>(m[i]);
        if (ref.hcount[i] != tgt.hcount[t]) {
            bad += 4;
        }
        if (ref.charges[i] && tgt.charges[t] && *ref.charges[i] != *tgt.charges[t]) {
            bad += 2;
        }
        for (const auto& [j, order] : ref.adj[i]) {
            if (j > static_cast<int>(i)) {
                const int other = tgt.bond(m[i], m[static_cast<std::size_t>(j)]).value_or(0);
                if (order != 0 && other != 0 && order != other) {
                    bad += 1;
                }
            }
        }
    }
    return bad;
}

std::optional<std::vector<int>> matchGraph(const MolGraph& ref, const MolGraph& tgt) {
    std::optional<std::pair<int, std::vector<int>>> best;
    Matcher(ref, tgt).run([&](const std::vector<int>& m) {
        std::pair<int, std::vector<int>> scored{matchScore(ref, tgt, m), m};
        if (!best || scored < *best) {
            best = std::move(scored);
        }
    });
    if (!best) {
        return std::nullopt;
    }
    return best->second;
}

}  // namespace mdbridge::chem
