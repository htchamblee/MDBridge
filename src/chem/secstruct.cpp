#include "mdbridge/chem/secstruct.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/text.hpp"
#include "mdbridge/io/pdb.hpp"

namespace mdbridge::chem {

namespace {

std::string pairText(const ResiduePair& p) {
    return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")";
}

// repr() of one character of the decoded text; a surrogate-escaped byte is passed as that byte,
// which pythonStringRepr shows as \udcXX
std::string charRepr(char32_t ch) {
    std::string text;
    if (ch >= 0xDC80 && ch <= 0xDCFF) {
        text.push_back(static_cast<char>(ch - 0xDC00));
    } else {
        core::appendUtf8(text, ch);
    }
    return core::pythonStringRepr(text);
}

// numpy's sort order for floats: NaN last
bool lessNanLast(double a, double b) {
    if (std::isnan(a)) {
        return false;
    }
    if (std::isnan(b)) {
        return true;
    }
    return a < b;
}

// np.linalg.norm(a - b, axis=-1): sqrt((dx*dx + dy*dy) + dz*dz); with `image`, of image(a - b)
double distance(const double* a, const double* b, const PairDetector::DifferenceImage* image = nullptr) {
    std::array<double, 3> v{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    if (image != nullptr) {
        v = (*image)(v);
    }
    const double dx = v[0];
    const double dy = v[1];
    const double dz = v[2];
    const double xx = dx * dx;
    const double yy = dy * dy;
    const double zz = dz * dz;
    const double s = xx + yy;
    return std::sqrt(s + zz);
}

}  // namespace

std::string_view kindLabel(int kind) {
    switch (kind) {
    case kWatsonCrick:
        return "Watson-Crick";
    case kWobble:
        return "wobble";
    case kOther:
        return "non-canonical";
    default:
        return "";
    }
}

const std::vector<PairRule>& pairRules() {
    static const std::vector<PairRule> rules = {
        {'G', "C", kWatsonCrick, {{"N1", "N3"}, {"N2", "O2"}, {"O6", "N4"}}, 2},
        {'A', "UT", kWatsonCrick, {{"N6", "O4"}, {"N1", "N3"}}, 2},
        {'G', "UT", kWobble, {{"O6", "N3"}, {"N1", "O2"}}, 2},
    };
    return rules;
}

const std::vector<std::string>& basePairAtoms() {
    static const std::vector<std::string> atoms = [] {
        std::set<std::string> names;
        for (const PairRule& rule : pairRules()) {
            for (const auto& [a, b] : rule.atoms) {
                names.emplace(a);
                names.emplace(b);
            }
        }
        return std::vector<std::string>(names.begin(), names.end());
    }();
    return atoms;
}

std::optional<char> baseLetter(std::string_view resname) {
    const auto info = io::parseResidueName(resname);
    if (!info) {
        return std::nullopt;
    }
    return info->base;
}

std::optional<std::string> canonicalName(std::string_view resname) {
    const auto info = io::parseResidueName(resname);
    if (!info) {
        return std::nullopt;
    }
    return info->sugar == "dna" ? std::string{'D', info->base} : std::string(1, info->base);
}

int lwKind(std::string_view code) {
    return code == "WCc" ? kWatsonCrick : code == "GUc" ? kWobble : kOther;
}

// --- dot-bracket

DotBracket dotbracketToPairs(std::string_view text) {
    const std::u32string chars = core::decodeUtf8Surrogateescape(text);
    std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> stacks(kOpeners.size());  // (n, pos)
    DotBracket out;
    std::int64_t n = 0;
    bool pendingBreak = false;
    std::int64_t pos = 0;
    for (const char32_t ch : chars) {
        ++pos;
        if (ch == static_cast<char32_t>(kChainSeparator) || core::isPythonSpace(ch)) {
            pendingBreak = n > 0;
            continue;
        }
        if (pendingBreak) {
            out.chainBreaks.push_back(n);
            pendingBreak = false;
        }
        const char c = ch < 0x80 ? static_cast<char>(ch) : '\0';
        const std::size_t opener = c ? kOpeners.find(c) : std::string_view::npos;
        const std::size_t closer = c ? kClosers.find(c) : std::string_view::npos;
        if (opener != std::string_view::npos) {
            stacks[opener].emplace_back(n, pos);
        } else if (closer != std::string_view::npos) {
            auto& stack = stacks[closer];
            if (stack.empty()) {
                throw core::ValueError("unbalanced dot-bracket: " + charRepr(ch) + " at position " +
                                       std::to_string(pos) + " (residue " + std::to_string(n + 1) +
                                       ") has no matching " + charRepr(static_cast<char32_t>(kOpeners[closer])));
            }
            out.pairs.emplace_back(stack.back().first, n);
            stack.pop_back();
        } else if (!c || kUnpaired.find(c) == std::string_view::npos) {
            throw core::ValueError("invalid dot-bracket character " + charRepr(ch) + " at position " +
                                   std::to_string(pos));
        }
        ++n;
    }
    std::size_t open = 0;
    std::int64_t firstPos = 0;
    std::size_t firstOpener = 0;
    for (std::size_t o = 0; o < stacks.size(); ++o) {
        for (const auto& entry : stacks[o]) {
            const std::int64_t p = entry.second;
            if (open == 0 || p < firstPos) {
                firstPos = p;
                firstOpener = o;
            }
            ++open;
        }
    }
    if (open > 0) {
        const std::string more = open > 1 ? " (and " + std::to_string(open - 1) + " more)" : "";
        throw core::ValueError("unbalanced dot-bracket: " +
                               charRepr(static_cast<char32_t>(kOpeners[firstOpener])) + " at position " +
                               std::to_string(firstPos) + " is never closed" + more);
    }
    std::sort(out.pairs.begin(), out.pairs.end());
    out.n = n;
    return out;
}

std::vector<std::vector<ResiduePair>> stems(std::span<const ResiduePair> pairs) {
    std::map<std::int64_t, std::int64_t> partner;
    for (const auto& [i, j] : pairs) {
        partner[i] = j;
    }
    std::vector<ResiduePair> sorted(pairs.begin(), pairs.end());
    std::sort(sorted.begin(), sorted.end());
    std::set<std::int64_t> seen;
    std::vector<std::vector<ResiduePair>> out;
    for (auto [i, j] : sorted) {
        if (seen.contains(i)) {
            continue;
        }
        std::vector<ResiduePair> stem;
        for (;;) {
            const auto it = partner.find(i);
            if (it == partner.end() || it->second != j || seen.contains(i)) {
                break;
            }
            stem.emplace_back(i, j);
            seen.insert(i);
            ++i;
            --j;
        }
        out.push_back(std::move(stem));
    }
    return out;
}

bool crosses(const ResiduePair& a, const ResiduePair& b) {
    return (a.first < b.first && b.first < a.second && a.second < b.second) ||
           (b.first < a.first && a.first < b.second && b.second < a.second);
}

std::string pairsToDotbracket(std::span<const ResiduePair> pairs, std::int64_t n,
                              std::span<const std::int64_t> chainBreaks) {
    std::vector<ResiduePair> norm;
    std::set<ResiduePair> inNorm;
    std::map<std::int64_t, ResiduePair> used;
    for (const ResiduePair& p : pairs) {
        const std::int64_t i = std::min(p.first, p.second);
        const std::int64_t j = std::max(p.first, p.second);
        if (i == j || i < 0 || j >= n) {
            throw core::ValueError("invalid base pair (" + std::to_string(p.first) + ", " +
                                   std::to_string(p.second) + ") for " + std::to_string(n) + " residues");
        }
        for (const std::int64_t k : {i, j}) {
            const auto it = used.find(k);
            if (it != used.end() && it->second != ResiduePair{i, j}) {
                throw core::ValueError("residue " + std::to_string(k) + " is in more than one base pair (" +
                                       pairText(it->second) + " and " + pairText({i, j}) + ")");
            }
            used[k] = {i, j};
        }
        if (inNorm.insert({i, j}).second) {
            norm.emplace_back(i, j);
        }
    }
    std::vector<std::vector<ResiduePair>> helices = stems(norm);
    // sorted by (-len(s), s[0]); every stem is non-empty here and starts at its own residue
    std::stable_sort(helices.begin(), helices.end(), [](const auto& a, const auto& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size();
        }
        return a.front() < b.front();
    });
    std::vector<std::vector<ResiduePair>> levels;
    std::string chars(static_cast<std::size_t>(std::max<std::int64_t>(n, 0)), '.');
    for (const auto& stem : helices) {
        std::size_t lev = 0;
        for (; lev < levels.size(); ++lev) {
            bool crossing = false;
            for (const ResiduePair& a : stem) {
                for (const ResiduePair& b : levels[lev]) {
                    if (crosses(a, b)) {
                        crossing = true;
                        break;
                    }
                }
                if (crossing) {
                    break;
                }
            }
            if (!crossing) {
                break;
            }
        }
        if (lev == levels.size()) {
            if (lev >= kOpeners.size()) {
                throw core::ValueError("too many pseudoknot levels for dot-bracket notation");
            }
            levels.emplace_back();
        }
        levels[lev].insert(levels[lev].end(), stem.begin(), stem.end());
        for (const auto& [i, j] : stem) {
            chars[static_cast<std::size_t>(i)] = kOpeners[lev];
            chars[static_cast<std::size_t>(j)] = kClosers[lev];
        }
    }
    std::set<std::int64_t> breaks;
    for (const std::int64_t b : chainBreaks) {
        if (0 < b && b < n) {
            breaks.insert(b);
        }
    }
    std::string out;
    out.reserve(chars.size() + breaks.size());
    for (std::size_t k = 0; k < chars.size(); ++k) {
        if (breaks.contains(static_cast<std::int64_t>(k))) {
            out += kChainSeparator;
        }
        out += chars[k];
    }
    return out;
}

UniquePairs uniquePairs(std::span<const ResiduePair> pairs) {
    std::set<std::int64_t> used;
    UniquePairs out;
    for (const ResiduePair& p : pairs) {
        if (used.contains(p.first) || used.contains(p.second)) {
            out.dropped.push_back(p);
            continue;
        }
        used.insert(p.first);
        used.insert(p.second);
        out.keep.push_back(p);
    }
    return out;
}

// --- builtin detection

namespace {

std::int64_t lookup(const AtomIndex& index, std::string_view name) {
    const auto it = index.find(name);
    return it == index.end() ? -1 : it->second;
}

std::int64_t c1Row(const AtomIndex& index) {
    const auto it = index.find("C1'");
    return it != index.end() ? it->second : lookup(index, "C1*");
}

}  // namespace

PairDetector::PairDetector(std::span<const std::string> resnames, std::span<const AtomIndex> atomIndex,
                           double cutoff, bool wobble)
    : m_n(resnames.size()), m_cutoff(cutoff) {
    if (atomIndex.size() < resnames.size()) {
        throw std::invalid_argument("PairDetector: one atom index per residue needed");
    }
    // MDATK 0.1.x: a CCD component's own atom names read as the standard ones (0DA's N1A is N1; an
    // atom of the standard name wins over a renamed one), and a name taken for a nucleotide only with
    // a backbone atom (the detergent LDA, a lanthanum ion LA are none)
    std::vector<AtomIndex> standardIndex(resnames.size());
    std::vector<std::optional<char>> bases;
    bases.reserve(resnames.size());
    for (std::size_t res = 0; res < resnames.size(); ++res) {
        AtomIndex& index = standardIndex[res];
        for (const auto& [name, row] : atomIndex[res]) {
            if (io::standardAtomName(resnames[res], name) == name) {
                index.emplace(name, row);
            }
        }
        for (const auto& [name, row] : atomIndex[res]) {
            index.emplace(io::standardAtomName(resnames[res], name), row);
        }
        std::vector<std::string_view> names;
        for (const auto& entry : index) {
            names.push_back(entry.first);
        }
        bases.push_back(io::hasBackboneAtoms(names) ? baseLetter(resnames[res]) : std::nullopt);
    }
    for (const PairRule& rule : pairRules()) {
        if (rule.kind == kWobble && !wobble) {
            continue;
        }
        Rule r;
        r.kind = rule.kind;
        r.required = rule.required;
        r.width = rule.atoms.size();
        const auto collect = [&](bool first, std::vector<std::int64_t>& rows, std::vector<std::int64_t>& atoms,
                                 std::vector<std::int64_t>& c1) {
            for (std::size_t res = 0; res < bases.size(); ++res) {
                const auto& b = bases[res];
                const bool match = b && (first ? *b == rule.base : rule.partners.find(*b) != std::string_view::npos);
                if (!match) {
                    continue;
                }
                std::vector<std::int64_t> idx;
                int present = 0;
                for (const auto& [pa, qa] : rule.atoms) {
                    idx.push_back(lookup(standardIndex[res], first ? pa : qa));
                    present += idx.back() >= 0 ? 1 : 0;
                }
                if (present < rule.required) {  // can no longer reach `required` bonds
                    continue;
                }
                rows.push_back(static_cast<std::int64_t>(res));
                atoms.insert(atoms.end(), idx.begin(), idx.end());
                c1.push_back(c1Row(standardIndex[res]));
            }
        };
        collect(true, r.rowsP, r.atomsP, r.c1P);
        collect(false, r.rowsQ, r.atomsQ, r.c1Q);
        if (!r.rowsP.empty() && !r.rowsQ.empty()) {
            m_rules.push_back(std::move(r));
        }
    }
}

std::vector<PairCandidate> PairDetector::candidates(std::span<const double> xyz, const DifferenceImage* image) const {
    if (xyz.size() % 3 != 0) {
        throw std::invalid_argument("PairDetector: coordinates must have 3 values per row");
    }
    const std::size_t nRows = xyz.size() / 3;
    const auto at = [&](std::int64_t row) -> const double* {
        const auto r = static_cast<std::size_t>(std::max<std::int64_t>(row, 0));
        if (r >= nRows) {
            throw std::out_of_range("PairDetector: atom row " + std::to_string(r) + " out of range for " +
                                    std::to_string(nRows) + " coordinates");
        }
        return xyz.data() + 3 * r;
    };
    std::vector<PairCandidate> out;
    std::vector<double> d;
    for (const Rule& rule : m_rules) {
        const std::size_t k = rule.width;
        // every row the rule reads must exist (numpy indexes them all)
        for (const auto* rows : {&rule.atomsP, &rule.atomsQ, &rule.c1P, &rule.c1Q}) {
            for (const std::int64_t row : *rows) {
                at(row);
            }
        }
        d.assign(k, 0.0);
        for (std::size_t a = 0; a < rule.rowsP.size(); ++a) {
            const std::int64_t* pa = rule.atomsP.data() + a * k;
            for (std::size_t b = 0; b < rule.rowsQ.size(); ++b) {
                const std::int64_t* qb = rule.atomsQ.data() + b * k;
                int bonds = 0;
                for (std::size_t t = 0; t < k; ++t) {
                    d[t] = (pa[t] < 0 || qb[t] < 0) ? HUGE_VAL : distance(at(pa[t]), at(qb[t]), image);
                    bonds += d[t] <= m_cutoff ? 1 : 0;
                }
                if (bonds < rule.required) {
                    continue;
                }
                const std::int64_t ca = rule.c1P[a];
                const std::int64_t cb = rule.c1Q[b];
                const bool both = ca >= 0 && cb >= 0;
                const double dc = distance(at(ca), at(cb), image);
                if (both && !(dc >= kC1WindowLow && dc <= kC1WindowHigh)) {
                    continue;
                }
                // np.sort(d)[:2].sum(): the two smallest (NaN last)
                std::vector<double> sorted = d;
                std::sort(sorted.begin(), sorted.end(), lessNanLast);
                PairCandidate c;
                c.first = rule.rowsP[a];
                c.second = rule.rowsQ[b];
                c.kind = rule.kind;
                c.score = sorted[0] + sorted[1];
                c.c1Deviation = both ? std::abs(dc - kC1Ideal) : kC1Ideal;
                out.push_back(c);
            }
        }
    }
    return out;
}

std::vector<BasePair> PairDetector::select(std::span<const PairCandidate> candidates, std::size_t nResidues) {
    const std::size_t m = candidates.size();
    std::vector<double> rounded(m);
    for (std::size_t k = 0; k < m; ++k) {
        // np.round(score / SCORE_TIE) * SCORE_TIE, round half to even
        rounded[k] = std::nearbyint(candidates[k].score / kScoreTie) * kScoreTie;
    }
    std::vector<std::size_t> order(m);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (lessNanLast(rounded[a], rounded[b])) {
            return true;
        }
        if (lessNanLast(rounded[b], rounded[a])) {
            return false;
        }
        return lessNanLast(candidates[a].c1Deviation, candidates[b].c1Deviation);
    });
    std::vector<bool> used(nResidues, false);
    std::vector<BasePair> keep;
    for (const std::size_t k : order) {
        const auto i = static_cast<std::size_t>(candidates[k].first);
        const auto j = static_cast<std::size_t>(candidates[k].second);
        if (i >= nResidues || j >= nResidues) {
            throw std::out_of_range("PairDetector: candidate residue out of range");
        }
        if (!used[i] && !used[j]) {
            used[i] = used[j] = true;
            keep.push_back({std::min(candidates[k].first, candidates[k].second),
                            std::max(candidates[k].first, candidates[k].second), candidates[k].kind});
        }
    }
    std::sort(keep.begin(), keep.end(), [](const BasePair& a, const BasePair& b) {
        return std::tie(a.i, a.j, a.kind) < std::tie(b.i, b.j, b.kind);
    });
    return keep;
}

std::vector<BasePair> PairDetector::detect(std::span<const double> xyz, const DifferenceImage* image) const {
    return select(candidates(xyz, image), m_n);
}

}  // namespace mdbridge::chem
