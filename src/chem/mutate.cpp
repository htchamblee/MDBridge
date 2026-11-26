// MDATK's base mutations (mdatk/builder.py v0.0.4.1 parse_mutations, mutated_name, _partners,
// _swap_base, mutate; MDATK's own code, ported) with CH-08 (MDATK 0.1.0's mutation specs) and CH-09
// (D11: ideal hydrogens on a mutated base). Arithmetic in NumPy's order (no fused multiply-adds:
// CMakeLists.txt); the Kabsch rotation (mdatk.spatial.kabsch) from core::svd3x3.

#include "mdbridge/chem/mutate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "mdbridge/chem/helix.hpp"
#include "mdbridge/chem/pdbtools.hpp"
#include "mdbridge/chem/secstruct.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/linalg.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

namespace {

using io::PdbAtom;
using io::PdbResidue;
using io::PdbStructure;
using Mat3 = std::array<Vec3, 3>;

constexpr std::string_view kReadHint = " (expected e.g. A:5G - chain, residue number, new base)";

// --- the base tables -----------------------------------------------------------------------------

struct Glycosidic {
    std::string_view n;
    std::string_view a;
    std::string_view b;
};
// GLYCOSIDIC
Glycosidic glycosidic(char base) {
    if (base == 'A' || base == 'G') {
        return {"N9", "C4", "C8"};
    }
    return {"N1", "C2", "C6"};
}

// _SUGAR_EXTRA
bool sugarExtra(std::string_view name) {
    static const std::set<std::string_view, std::less<>> names = {"P",    "OP1",  "OP2",  "OP3", "O1P", "O2P",
                                                                  "O3P",  "HOP2", "HOP3", "H5T", "H3T"};
    return names.count(name) > 0;
}

// _ALL_BASE_NAMES: every base atom and base hydrogen of A C G U, plus the names other programs give
bool anyBaseName(std::string_view name) {
    static const std::set<std::string, std::less<>> names = [] {
        std::set<std::string, std::less<>> out = {"H62", "H22", "H42", "C7", "H71", "H72", "H73"};
        for (const char b : std::string_view("ACGU")) {
            for (const std::string_view n : baseAtomNames(b)) {
                out.emplace(n);
            }
            for (const std::string_view n : baseHydrogenNames(b)) {
                out.emplace(n);
            }
        }
        return out;
    }();
    return names.count(name) > 0;
}

bool contains(std::span<const std::string_view> names, std::string_view name) {
    return std::find(names.begin(), names.end(), name) != names.end();
}

// WATSON_CRICK_PAIRS, WOBBLE_PAIRS, COMPLEMENT
bool watsonCrickPair(char a, char b) {
    const std::string p{a, b};
    return p == "AU" || p == "UA" || p == "GC" || p == "CG";
}
bool wobblePair(char a, char b) {
    const std::string p{a, b};
    return p == "GU" || p == "UG";
}
char complementOf(char b) {
    switch (b) {
        case 'A':
            return 'U';
        case 'U':
            return 'A';
        case 'G':
            return 'C';
        default:
            return 'G';
    }
}

// _new_partner_base
char newPartnerBase(char base, char current, bool wobble) {
    if (wobble && wobblePair(base, current)) {
        return current;
    }
    return complementOf(base);
}

// _element: the name without leading digits, its first letter
std::string elementOf(std::string_view name) {
    const std::size_t k = name.find_first_not_of("0123456789");
    return k == std::string_view::npos ? std::string() : std::string(name.substr(k, 1));
}

std::string starToPrime(std::string name) {
    std::replace(name.begin(), name.end(), '*', '\'');
    return name;
}

const Vec3* findCoord(const AtomCoords& coords, std::string_view name) {
    for (const auto& [n, xyz] : coords) {
        if (n == name) {
            return &xyz;
        }
    }
    return nullptr;
}

std::string fixed2(double x) {
    char text[64];
    std::snprintf(text, sizeof text, "%.2f", x);
    return text;
}

std::string join(const std::vector<std::string>& items, std::string_view sep) {
    std::string out;
    for (std::size_t k = 0; k < items.size(); ++k) {
        if (k) {
            out += sep;
        }
        out += items[k];
    }
    return out;
}

std::string optionalCall(const std::optional<char>& call) {
    return call ? std::string(1, *call) : std::string("None");
}

// --- builder._nucleotide -------------------------------------------------------------------------

// parse_residue_name of the name, guarded by a backbone atom (a lanthanum ion LA is no RNA)
std::optional<io::ResidueName> builderNucleotide(const PdbResidue& res) {
    std::optional<io::ResidueName> info = io::parseResidueName(res.name);
    if (!info) {
        return std::nullopt;
    }
    for (const PdbAtom& a : res.atoms) {
        if (io::isBackboneAtomName(a.name)) {
            return info;
        }
    }
    return std::nullopt;
}

// --- kabsch (mdatk.spatial.kabsch with weights) -----------------------------------------------

struct Kabsch {
    Mat3 rotation{};
    Vec3 mobileCentre{};
    Vec3 referenceCentre{};
};

Kabsch kabsch(std::span<const Vec3> mobile, std::span<const Vec3> reference, std::span<const double> weights) {
    const std::size_t n = mobile.size();
    double total = 0.0;
    for (const double weight : weights) {
        total += weight;
    }
    std::vector<double> w(n);
    for (std::size_t k = 0; k < n; ++k) {
        w[k] = weights[k] / total;
    }
    Kabsch out;
    out.rotation = {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    for (std::size_t k = 0; k < n; ++k) {  // w @ mobile, w @ reference
        for (std::size_t c = 0; c < 3; ++c) {
            out.mobileCentre[c] += w[k] * mobile[k][c];
            out.referenceCentre[c] += w[k] * reference[k][c];
        }
    }
    if (n < 3) {
        return out;
    }
    // h = ((mobile - mc) * w[:, None]).T @ (reference - rc)
    core::Mat3 h{};
    for (std::size_t k = 0; k < n; ++k) {
        Vec3 x{};
        Vec3 y{};
        for (std::size_t c = 0; c < 3; ++c) {
            x[c] = (mobile[k][c] - out.mobileCentre[c]) * w[k];
            y[c] = reference[k][c] - out.referenceCentre[c];
        }
        for (std::size_t p = 0; p < 3; ++p) {
            for (std::size_t q = 0; q < 3; ++q) {
                h[p][q] += x[p] * y[q];
            }
        }
    }
    for (const auto& row : h) {
        for (const double v : row) {
            if (!std::isfinite(v)) {
                throw BuildError("the glycosidic-frame fit failed (non-finite coordinates)");
            }
        }
    }
    const core::Svd3 svd = core::svd3x3(h);
    core::Mat3 vu{};  // vt.T @ u.T
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            vu[i][j] = svd.vt[0][i] * svd.u[j][0] + svd.vt[1][i] * svd.u[j][1] + svd.vt[2][i] * svd.u[j][2];
        }
    }
    const double det = core::det3(vu);
    const double d = det > 0 ? 1.0 : det < 0 ? -1.0 : 1.0;  // np.sign, 0 -> 1
    for (std::size_t i = 0; i < 3; ++i) {  // vt.T @ diag(1, 1, d) @ u.T
        const Vec3 v{svd.vt[0][i], svd.vt[1][i], svd.vt[2][i] * d};
        for (std::size_t j = 0; j < 3; ++j) {
            out.rotation[i][j] = v[0] * svd.u[j][0] + v[1] * svd.u[j][1] + v[2] * svd.u[j][2];
        }
    }
    return out;
}

// (x - mc) @ rot.T + rc
Vec3 moved(const Kabsch& k, const Vec3& x) {
    const Vec3 v{x[0] - k.mobileCentre[0], x[1] - k.mobileCentre[1], x[2] - k.mobileCentre[2]};
    Vec3 out{};
    for (std::size_t j = 0; j < 3; ++j) {
        double sum = v[0] * k.rotation[j][0];
        sum += v[1] * k.rotation[j][1];
        sum += v[2] * k.rotation[j][2];
        out[j] = sum + k.referenceCentre[j];
    }
    return out;
}

// --- mutation specs ------------------------------------------------------------------------------

bool asciiLetter(char32_t c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}
bool baseLetter(char32_t c) {
    return std::u32string_view(U"ACGUacgu").find(c) != std::u32string_view::npos;
}

struct SpecMatch {
    std::u32string chain;
    std::u32string number;
    std::u32string icode;
    char32_t base = 0;
};

// (-?\d+)([A-Za-z]?)(?:>|:|->)?([ACGUacgu])$ from `p`
std::optional<SpecMatch> matchRest(const std::u32string& s, std::size_t p) {
    SpecMatch m;
    std::size_t q = p;
    if (q < s.size() && s[q] == U'-') {
        ++q;
    }
    const std::size_t digits = q;
    while (q < s.size() && core::isPythonDigit(s[q])) {
        ++q;
    }
    if (q == digits) {
        return std::nullopt;
    }
    m.number = s.substr(p, q - p);
    const auto tail = [&](std::size_t r) -> std::optional<char32_t> {  // separator? base $
        for (const std::u32string_view sep : {std::u32string_view(U">"), std::u32string_view(U":"),
                                              std::u32string_view(U"->"), std::u32string_view()}) {
            if (s.compare(r, sep.size(), sep) != 0 && !sep.empty()) {
                continue;
            }
            const std::size_t b = r + sep.size();
            if (b + 1 == s.size() && baseLetter(s[b])) {
                return s[b];
            }
        }
        return std::nullopt;
    };
    if (q < s.size() && asciiLetter(s[q])) {
        if (const std::optional<char32_t> base = tail(q + 1)) {
            m.icode = s.substr(q, 1);
            m.base = *base;
            return m;
        }
    }
    if (const std::optional<char32_t> base = tail(q)) {
        m.base = *base;
        return m;
    }
    return std::nullopt;
}

// _MUTATION_RE of MDATK 0.1.0 (A C G U as in 0.0.4.1):
// ^(?:([^:\s]?):|([^:\s\d\-]?))(-?\d+)([A-Za-z]?)(?:>|:|->)?([ACGUacgu])$
std::optional<SpecMatch> matchSpec(const std::u32string& s) {
    const auto chainChar = [](char32_t c) { return c != U':' && !core::isPythonSpace(c); };
    std::vector<std::pair<std::u32string, std::size_t>> prefixes;  // (chain, rest), in the regex's order
    if (s.size() >= 2 && chainChar(s[0]) && s[1] == U':') {
        prefixes.emplace_back(s.substr(0, 1), 2);
    }
    if (!s.empty() && s[0] == U':') {
        prefixes.emplace_back(std::u32string(), 1);
    }
    if (!s.empty() && chainChar(s[0]) && !core::isPythonDigit(s[0]) && s[0] != U'-') {
        prefixes.emplace_back(s.substr(0, 1), 1);
    }
    prefixes.emplace_back(std::u32string(), 0);
    for (const auto& [chain, rest] : prefixes) {
        if (std::optional<SpecMatch> m = matchRest(s, rest)) {
            m->chain = chain;
            return m;
        }
    }
    return std::nullopt;
}

void addMutation(std::vector<Mutation>& out, const std::string& item) {
    const std::string stripped = core::pythonStrip(item);
    const std::optional<SpecMatch> m = matchSpec(core::decodeUtf8(stripped));
    const auto unreadable = [&] {
        return BuildError("cannot read mutation " + core::pythonStringRepr(item) + std::string(kReadHint));
    };
    if (!m) {
        throw unreadable();
    }
    Mutation mut;
    mut.chain = core::encodeUtf8(m->chain);
    try {
        mut.resseq = core::pythonInt(core::encodeUtf8(m->number));
    } catch (const std::overflow_error&) {
        throw unreadable();  // beyond 64 bits: no PDB residue has such a number
    }
    mut.icode = core::encodeUtf8(m->icode);
    mut.base = static_cast<char>(m->base >= U'a' ? m->base - U'a' + U'A' : m->base);
    for (const Mutation& o : out) {
        if (o.chain == mut.chain && o.resseq == mut.resseq && o.icode == mut.icode) {
            throw BuildError("residue " + mut.where() + " is mutated twice");
        }
    }
    out.push_back(std::move(mut));
}

// re.split(r"[,;\s]+", spec.strip())
std::vector<std::string> splitItems(std::string_view spec) {
    const std::u32string s = core::decodeUtf8(core::pythonStrip(spec));
    std::vector<std::string> out;
    std::u32string current;
    bool inSeparator = false;
    for (const char32_t c : s) {
        if (c == U',' || c == U';' || core::isPythonSpace(c)) {
            if (!inSeparator) {
                out.push_back(core::encodeUtf8(current));
                current.clear();
            }
            inSeparator = true;
        } else {
            current.push_back(c);
            inSeparator = false;
        }
    }
    out.push_back(core::encodeUtf8(current));
    return out;
}

std::vector<Mutation> parseItems(const std::vector<std::string>& items) {
    std::vector<Mutation> out;
    for (const std::string& item : items) {
        if (!item.empty()) {
            addMutation(out, item);
        }
    }
    if (out.empty()) {
        throw BuildError("no mutations given (e.g. A:5G, B:20C)");
    }
    return out;
}

// --- base pairs in the input (_partners) -------------------------------------------------------

// residue index -> index of its partner (canonical pairs of the builtin detector, G·U included)
std::map<std::size_t, std::size_t> partners(const std::vector<PdbResidue*>& residues, double cutoff) {
    std::vector<std::string> names;
    std::vector<AtomIndex> index;
    std::vector<double> xyz;
    std::int64_t n = 0;
    for (const PdbResidue* r : residues) {
        AtomIndex rows;
        for (const auto& [name, p] : r->coords()) {
            rows.emplace(name, n++);
            xyz.insert(xyz.end(), p.begin(), p.end());
        }
        index.push_back(std::move(rows));
        names.push_back(builderNucleotide(*r) ? r->name : std::string());
    }
    const PairDetector detector(names, index, cutoff, true);
    std::map<std::size_t, std::size_t> out;
    for (const BasePair& p : detector.detect(xyz)) {
        out[static_cast<std::size_t>(p.i)] = static_cast<std::size_t>(p.j);
        out[static_cast<std::size_t>(p.j)] = static_cast<std::size_t>(p.i);
    }
    return out;
}

// --- _swap_base ------------------------------------------------------------------------------------

std::string whereOf(const PdbResidue& res) {
    return (res.chain.empty() ? std::string("_") : res.chain) + ":" + std::to_string(res.resseq) + res.icode;
}

core::Value volumesValue(const SugarChirality& c) {
    core::Table t;
    for (const auto& [name, v] : c.volumes) {
        t.emplace_back(name, core::Value(v));
    }
    return core::Value(std::move(t));
}

// Replace the base of `res` in place; returns the notes
Notes swapBase(PdbResidue& res, char base) {
    const char old = io::parseResidueName(res.name)->base;
    const AtomCoords coords = res.coords();
    const Glycosidic go = glycosidic(old);
    const std::array<std::string_view, 4> oldFrame = {"C1'", go.n, go.a, go.b};
    std::vector<std::string> missing;
    for (const std::string_view a : oldFrame) {
        if (!findCoord(coords, a)) {
            missing.emplace_back(a);
        }
    }
    if (!missing.empty()) {
        throw BuildError(res.label() + ": missing atom(s) " + join(missing, ", ") + " of the glycosidic frame");
    }
    const std::optional<char> before = sugarChirality(coords, old).call;
    const std::vector<TemplateAtom> t = helixTemplate(std::string_view(&base, 1), true);  // CH-09: with hydrogens
    const auto templ = [&t](std::string_view name) -> const Vec3& {
        for (const TemplateAtom& a : t) {
            if (a.name == name) {
                return a.xyz;
            }
        }
        throw std::logic_error("mutate: no template atom " + std::string(name));
    };
    const Glycosidic gn = glycosidic(base);
    const std::array<std::string_view, 4> newFrame = {"C1'", gn.n, gn.a, gn.b};
    std::array<Vec3, 4> mob{};
    std::array<Vec3, 4> ref{};
    for (std::size_t k = 0; k < 4; ++k) {
        mob[k] = templ(newFrame[k]);
        ref[k] = *findCoord(coords, oldFrame[k]);
    }
    const std::array<double, 4> weights = {2.0, 2.0, 1.0, 1.0};  // the glycosidic bond itself counts double
    const Kabsch fit = kabsch(mob, ref, weights);
    double sum = 0.0;  // ((moved - ref) ** 2).sum(axis=1).mean()
    for (std::size_t k = 0; k < 4; ++k) {
        const Vec3 m = moved(fit, mob[k]);
        double row = (m[0] - ref[k][0]) * (m[0] - ref[k][0]);
        row += (m[1] - ref[k][1]) * (m[1] - ref[k][1]);
        row += (m[2] - ref[k][2]) * (m[2] - ref[k][2]);
        sum += row;
    }
    const double rms = std::sqrt(sum / 4.0);
    // record, B-factor and segid of the new base atoms come from the C1' the fit used (coords() takes
    // the first alternate location); the new base is one conformer with full occupancy
    const PdbAtom* c1 = &res.atoms.front();
    for (const PdbAtom& a : res.atoms) {
        if (starToPrime(a.name) == "C1'") {
            c1 = &a;
            break;
        }
    }
    const PdbAtom templateC1 = *c1;
    std::set<std::string> altSet;
    for (const PdbAtom& a : res.atoms) {
        if (!a.alt.empty()) {
            altSet.insert(a.alt);
        }
    }
    const std::vector<std::string> altlocs(altSet.begin(), altSet.end());  // sorted
    std::vector<PdbAtom> keep;
    std::vector<std::string> droppedH;
    std::vector<std::string> unknown;
    std::optional<std::size_t> insertAt;
    for (const PdbAtom& a : res.atoms) {
        const std::string name = starToPrime(a.name);
        const bool isSugar = name.find('\'') != std::string::npos || sugarExtra(name);
        if (isSugar) {
            keep.push_back(a);
        } else if (anyBaseName(name) || contains(baseAtomNames(old), name)) {
            if (!insertAt) {
                insertAt = keep.size();
            }
            if (elementOf(name) == "H") {
                droppedH.push_back(name);
            }
        } else {
            unknown.push_back(name);
            keep.push_back(a);
        }
    }
    const std::size_t at = insertAt.value_or(keep.size());
    std::vector<std::string> newNames;
    for (const std::string_view n : baseAtomNames(base)) {
        newNames.emplace_back(n);
    }
    std::vector<std::string> addedH;  // CH-09
    if (!droppedH.empty()) {
        for (const std::string_view h : baseHydrogenNames(base)) {
            addedH.emplace_back(h);
        }
    }
    std::vector<PdbAtom> newAtoms;
    for (const std::vector<std::string>* names : {&newNames, &addedH}) {
        for (const std::string& name : *names) {
            const Vec3 xyz = moved(fit, templ(name));
            PdbAtom atom;
            atom.record = templateC1.record;
            atom.name = name;
            // occupancy 1.00, the C1' B-factor and the rest of its columns
            atom.tail = templateC1.tail.substr(0, 3) + io::newAtomTail(xyz).substr(3, 30) + templateC1.tail.substr(33);
            atom.element = elementOf(name);
            newAtoms.push_back(std::move(atom));
        }
    }
    std::vector<PdbAtom> atoms(keep.begin(), keep.begin() + static_cast<std::ptrdiff_t>(at));
    atoms.insert(atoms.end(), newAtoms.begin(), newAtoms.end());
    atoms.insert(atoms.end(), keep.begin() + static_cast<std::ptrdiff_t>(at), keep.end());
    res.atoms = std::move(atoms);
    const std::string oldName = res.name;
    res.name = mutatedName(res.name, base);
    const SugarChirality after = sugarChirality(res.coords(), base);
    if (before && after.call != before) {
        throw BuildError(res.label() + ": the base swap changed the sugar chirality call (" + optionalCall(before) +
                         " -> " + optionalCall(after.call) + ", " + core::pythonRepr(volumesValue(after)) + ")");
    }
    const std::string where = whereOf(res);
    Notes notes = {where + " " + std::string(1, old) + "->" + std::string(1, base) + " (" + oldName + " -> " +
                   res.name + ", glycosidic-frame fit " + fixed2(rms) + " Å)"};
    if (!droppedH.empty()) {
        notes.push_back(where + ": replaced " + std::to_string(droppedH.size()) + " hydrogen(s) of the old base (" +
                        join(droppedH, " ") + ") with " + std::to_string(addedH.size()) + " of the new base (" +
                        join(addedH, " ") + ") at ideal positions");
    }
    if (!unknown.empty()) {
        notes.push_back("WARNING: " + res.label() + ": kept unrecognised atom(s) " + join(unknown, " "));
    }
    if (!altlocs.empty()) {
        notes.push_back(where + ": alternate locations " + join(altlocs, ",") + "; the new base is fitted on " +
                        (templateC1.alt.empty() ? altlocs.front() : templateC1.alt) +
                        " and written as one conformer with occupancy 1.00");
    }
    return notes;
}

MutateResult mutateStructure(PdbStructure s, const std::vector<Mutation>& mutations, const MutateOptions& options) {
    Notes notes = s.notes;
    const Notes chirNotes = assignChirality(s, true);
    std::vector<PdbResidue*> residues;
    for (std::vector<PdbResidue>& seg : s.segments) {
        for (PdbResidue& r : seg) {
            residues.push_back(&r);
        }
    }
    std::map<std::tuple<std::string, std::int64_t, std::string>, std::size_t> byKey;  // the last one wins
    for (std::size_t k = 0; k < residues.size(); ++k) {
        byKey[{residues[k]->chain, residues[k]->resseq, residues[k]->icode}] = k;
    }
    std::vector<std::pair<std::size_t, Mutation>> targets;
    for (const Mutation& m : mutations) {
        const auto hit = byKey.find({m.chain, m.resseq, m.icode});
        if (hit == byKey.end()) {
            throw BuildError("residue " + m.where() + " not found");
        }
        const std::optional<io::ResidueName> info = builderNucleotide(*residues[hit->second]);
        if (!info || info->sugar != "rna") {
            throw BuildError(residues[hit->second]->label() + " is not an RNA nucleotide");
        }
        targets.emplace_back(hit->second, m);
    }
    if (options.compensatory) {
        const std::map<std::size_t, std::size_t> partnerOf = partners(residues, options.cutoff);
        std::set<std::size_t> explicitTargets;
        for (const auto& target : targets) {
            explicitTargets.insert(target.first);
        }
        const std::size_t given = targets.size();
        for (std::size_t t = 0; t < given; ++t) {
            const std::size_t k = targets[t].first;
            const Mutation m = targets[t].second;
            const auto hit = partnerOf.find(k);
            if (hit == partnerOf.end()) {
                notes.push_back(m.where() + ": no Watson-Crick / G·U partner in the input; no compensatory change");
                continue;
            }
            const PdbResidue& other = *residues[hit->second];
            const std::string otherAt = other.chain + ":" + std::to_string(other.resseq);
            const std::optional<io::ResidueName> otherInfo = builderNucleotide(other);
            if (!otherInfo || otherInfo->sugar != "rna") {
                notes.push_back(m.where() + ": its partner " + otherAt + " " + other.name +
                                " is not an RNA nucleotide; no compensatory change");
                continue;
            }
            const char otherBase = otherInfo->base;
            if (explicitTargets.count(hit->second)) {
                char newOther = 0;
                for (const auto& [kk, mm] : targets) {
                    if (kk == hit->second) {
                        newOther = mm.base;
                        break;
                    }
                }
                if (!watsonCrickPair(m.base, newOther) && !wobblePair(m.base, newOther)) {
                    notes.push_back("WARNING: " + m.where() + " " + std::string(1, m.base) + " and its partner " +
                                    otherAt + " " + std::string(1, newOther) +
                                    " (both mutated explicitly) do not pair canonically");
                }
                continue;
            }
            const char newBase = newPartnerBase(m.base, otherBase, options.wobble);
            if (newBase == otherBase) {
                notes.push_back(m.where() + " " + std::string(1, m.base) + ": partner " + otherAt + " " +
                                std::string(1, otherBase) + " already pairs (" +
                                (wobblePair(m.base, otherBase) ? "G·U wobble" : "Watson-Crick") + ")");
                continue;
            }
            Mutation comp;
            comp.chain = other.chain;
            comp.resseq = other.resseq;
            comp.icode = other.icode;
            comp.base = newBase;
            comp.isExplicit = false;
            comp.reason = "compensates " + m.where() + " " + std::string(1, m.base);
            targets.emplace_back(hit->second, std::move(comp));
            explicitTargets.insert(hit->second);
        }
    }
    std::size_t done = 0;
    for (const auto& [k, m] : targets) {
        PdbResidue& res = *residues[k];
        if (io::parseResidueName(res.name)->base == m.base) {
            notes.push_back(m.where() + ": already " + std::string(1, m.base) + "; unchanged");
            continue;
        }
        Notes change = swapBase(res, m.base);
        if (!m.isExplicit) {
            change[0] += " - compensatory, " + m.reason;
        }
        notes.insert(notes.end(), change.begin(), change.end());
        ++done;
    }
    std::size_t stale = 0;
    std::vector<std::string> header;
    for (const std::string& line : s.header) {
        if (line.rfind("REMARK 950", 0) == 0) {
            ++stale;
        } else {
            header.push_back(line);
        }
    }
    if (stale) {
        s.header = std::move(header);
        notes.push_back("Dropped " + std::to_string(stale) +
                        " REMARK 950 base-pair annotation lines (they describe the old bases)");
    }
    if (!s.conect.empty()) {
        notes.push_back("Stripped " + std::to_string(s.conect.size()) + " CONECT records");
    }
    notes.push_back("Mutated " + std::to_string(done) + " residue(s)");
    for (const std::string& n : chirNotes) {
        if (n.rfind("WARNING", 0) == 0) {
            notes.push_back(n);
        }
    }
    return {chem::writePdb(s), std::move(notes)};
}

}  // namespace

// --- public ----------------------------------------------------------------------------------------

std::string Mutation::where() const {
    return (chain.empty() ? std::string("_") : chain) + ":" + std::to_string(resseq) + icode;
}

std::vector<Mutation> parseMutations(std::string_view spec) {
    return parseItems(splitItems(spec));
}

std::vector<Mutation> parseMutations(std::span<const std::string> items) {
    return parseItems(std::vector<std::string>(items.begin(), items.end()));
}

std::string mutatedName(std::string_view name, char base) {
    const std::string up = io::pythonUpper(core::pythonStrip(name));
    static const std::map<char, std::string_view> aliases = {{'A', "ADE"}, {'C', "CYT"}, {'G', "GUA"}, {'U', "URA"}};
    if (up == "ADE" || up == "CYT" || up == "GUA" || up == "URA" || up == "URI") {
        return std::string(aliases.at(base));
    }
    // ^([RL0]?)([ACGU])([53N]?)$
    std::size_t k = 0;
    std::string prefix;
    if (k < up.size() && (up[k] == 'R' || up[k] == 'L' || up[k] == '0') && up.size() >= 2 &&
        std::string_view("ACGU").find(up[1]) != std::string_view::npos) {
        prefix = up.substr(0, 1);
        k = 1;
    }
    if (k < up.size() && std::string_view("ACGU").find(up[k]) != std::string_view::npos) {
        const std::string rest = up.substr(k + 1);
        if (rest.empty() || rest == "5" || rest == "3" || rest == "N") {
            return prefix + base + rest;
        }
    }
    throw BuildError("cannot rename residue " + core::pythonStringRepr(std::string(name)) +
                     " (RNA nucleotides only)");
}

MutateResult mutate(const PdbSource& source, std::string_view mutations, const MutateOptions& options) {
    // MDATK reads the structure before the specs: a file that cannot be read fails first
    PdbStructure s = source.read();
    return mutateStructure(std::move(s), parseMutations(mutations), options);
}

MutateResult mutate(const PdbSource& source, std::vector<Mutation> mutations, const MutateOptions& options) {
    return mutateStructure(source.read(), mutations, options);
}


}  // namespace mdbridge::chem
