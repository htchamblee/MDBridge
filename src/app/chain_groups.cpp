#include "chain_groups.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string_view>

#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/io/loader.hpp"
#include "mdbridge/io/parm7.hpp"

namespace mdbridge::app {

namespace {

// residue names of water models (also recognised by composition: one oxygen, hydrogens, virtual sites)
const std::set<std::string, std::less<>> kWaterResnames{
    "WAT", "HOH", "H2O", "DOD", "SOL", "TIP", "TIP3", "TIP4", "TIP5", "TP3", "TP4", "TP5", "T3P", "T4P",
    "T5P", "T4E", "SPC", "SPCE", "OPC", "OPC3", "PL3", "WTR", "TIP3P", "TIP4P", "TIP5P"};
// MDATK's ion residue names (selections._ION_ELEMENTS) and a few more
const std::set<std::string, std::less<>> kIonResnames{
    "MG", "MG2", "NA", "NA+", "SOD", "K", "K+", "POT", "CL", "CL-", "CLA", "CA", "CAL", "ZN", "MN", "FE", "LI",
    "RB", "CS", "BR", "MN2", "ZN2", "CA2", "CD", "CO", "NI", "CU", "SR", "BA", "IOD"};
constexpr double kLinkMax = 2.0;  // largest heavy-atom distance (Å) that bonds two consecutive residues
// distinct colours, first the default groups' blue and orange; then golden-angle hues
constexpr const char* kChainColors[] = {"#2a78d6", "#eb6834", "#2f9e6e", "#c2408f", "#7a55d3",
                                        "#18a0b4", "#8c5a2b", "#5e7391", "#a3b02a", "#d4789c"};

std::string strip(std::string_view text) {
    const auto space = [](char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v'; };
    std::size_t a = 0;
    std::size_t b = text.size();
    while (a < b && space(text[a])) ++a;
    while (b > a && space(text[b - 1])) --b;
    return std::string(text.substr(a, b - a));
}

std::string upper(std::string text) {
    for (char& c : text) {
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    }
    return text;
}

bool isAlpha(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }
bool isDigit(char c) { return c >= '0' && c <= '9'; }

// [A-Za-z][A-Za-z0-9]*: the chain qualifiers residue ranges accept
bool rangeChain(std::string_view s) {
    if (s.empty() || !isAlpha(s[0])) return false;
    return std::all_of(s.begin(), s.end(), [](char c) { return isAlpha(c) || isDigit(c); });
}

// [A-Za-z0-9_]+
bool selectionWord(std::string_view s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return isAlpha(c) || isDigit(c) || c == '_'; });
}

// re.sub(r"[^A-Za-z0-9_.\-]", "_", text): one '_' per character (a UTF-8 sequence is one character)
std::string nameSafe(std::string_view text) {
    std::string out;
    for (std::size_t k = 0; k < text.size(); ++k) {
        const char c = text[k];
        if (isAlpha(c) || isDigit(c) || c == '_' || c == '.' || c == '-') {
            out += c;
        } else if ((static_cast<unsigned char>(c) & 0xC0) != 0x80) {
            out += '_';
        }
    }
    return out;
}

bool usable(const std::vector<std::string>& labels) {
    if (labels.empty() || !std::all_of(labels.begin(), labels.end(), [](const std::string& s) { return !s.empty(); })) {
        return false;
    }
    const std::set<std::string> distinct(labels.begin(), labels.end());
    return distinct.size() > 1 && !distinct.count("SYSTEM");
}

std::string freeName(const std::string& base, const std::set<std::string>& taken) {
    std::string name = base;
    for (int k = 2; taken.count(name); ++k) {
        name = base + "_" + std::to_string(k);
    }
    return name;
}

std::vector<std::pair<std::int64_t, std::int64_t>> runs(std::vector<std::int64_t> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    std::vector<std::pair<std::int64_t, std::int64_t>> out;
    for (const std::int64_t v : values) {
        if (!out.empty() && v == out.back().second + 1) {
            out.back().second = v;
        } else {
            out.emplace_back(v, v);
        }
    }
    return out;
}

std::string joinRuns(const std::vector<std::pair<std::int64_t, std::int64_t>>& rs, const char* dash,
                     const char* separator, const std::string& prefix = {}) {
    std::string out;
    for (std::size_t k = 0; k < rs.size(); ++k) {
        if (k) out += separator;
        out += prefix + std::to_string(rs[k].first);
        if (rs[k].first != rs[k].second) {
            out += dash + std::to_string(rs[k].second);
        }
    }
    return out;
}

// colorsys.hsv_to_rgb
std::array<double, 3> hsvToRgb(double h, double s, double v) {
    if (s == 0.0) return {v, v, v};
    int i = static_cast<int>(h * 6.0);
    const double f = h * 6.0 - i;
    const double p = v * (1.0 - s);
    const double q = v * (1.0 - s * f);
    const double t = v * (1.0 - s * (1.0 - f));
    i %= 6;
    switch (i) {
    case 0: return {v, t, p};
    case 1: return {q, v, p};
    case 2: return {p, v, t};
    case 3: return {p, q, v};
    case 4: return {t, p, v};
    default: return {v, p, q};
    }
}

// Colours not in `used` (case-insensitive): the palette first, then golden-angle hues.
class Colours {
public:
    explicit Colours(const std::vector<core::ResidueGroup>& used) {
        for (const auto& g : used) m_taken.insert(upper(strip(g.color)));
    }
    std::string next() {
        while (m_palette < std::size(kChainColors)) {
            const std::string c = kChainColors[m_palette++];
            if (!m_taken.count(upper(c))) return c;
        }
        for (;;) {
            const double h = std::fmod(0.11 + 0.618034 * static_cast<double>(m_hue++), 1.0);
            const auto rgb = hsvToRgb(h, 0.62, 0.78);
            char text[16];
            std::snprintf(text, sizeof text, "#%02x%02x%02x", static_cast<int>(rgb[0] * 255),
                          static_cast<int>(rgb[1] * 255), static_cast<int>(rgb[2] * 255));
            if (!m_taken.count(upper(text))) return text;
        }
    }

private:
    std::set<std::string> m_taken;
    std::size_t m_palette = 0;
    std::int64_t m_hue = 0;
};

std::optional<std::vector<std::int32_t>> atomsOf(const std::shared_ptr<const core::Topology>& topology,
                                                 const core::Frame* current, const core::ResidueGroup& group,
                                                 const std::vector<core::ResidueGroup>& groups) {
    try {
        return core::select(topology, current, "group:" + group.name, groups, true).indices();
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::vector<std::int32_t> atomsOfResidues(const core::Topology& top, const std::vector<std::int32_t>& residues) {
    std::vector<std::int32_t> out;
    for (const std::int32_t r : residues) {
        const auto atoms = top.residueAtoms(r);
        out.insert(out.end(), atoms.begin(), atoms.end());
    }
    std::sort(out.begin(), out.end());
    return out;
}

// (residues, selection) specs for a chain, simplest first
std::vector<std::pair<std::string, std::string>> specCandidates(const core::Topology& top, const PolymerChain& chain) {
    std::vector<std::int64_t> resids;
    for (const std::int32_t r : chain.residues) resids.push_back(top.resids[static_cast<std::size_t>(r)]);
    const auto rs = runs(resids);
    std::vector<std::pair<std::string, std::string>> out;
    out.emplace_back(joinRuns(rs, "-", ", "), "");
    const char* qualifierReads = top.chainIDs ? "chainID" : "segid";  // see rangesToSelection
    if (chain.source == qualifierReads && rangeChain(chain.label)) {
        out.emplace_back(joinRuns(rs, "-", ", ", chain.label + ":"), "");
    }
    if ((chain.source == "chainID" || chain.source == "segid") && selectionWord(chain.label)) {
        out.emplace_back("", chain.source + " " + chain.label + " and resid " + joinRuns(rs, ":", " "));
    }
    std::vector<std::int64_t> atoms;
    for (const std::int32_t a : atomsOfResidues(top, chain.residues)) atoms.push_back(a);
    out.emplace_back("", "index " + joinRuns(runs(atoms), ":", " "));
    return out;
}

}  // namespace

std::string ChainMerge::summary(const std::string& source) const {
    std::vector<std::string> parts;
    const auto join = [](const std::vector<std::string>& items) {
        std::string out;
        for (std::size_t k = 0; k < items.size(); ++k) out += (k ? ", " : "") + items[k];
        return out;
    };
    if (!replaced.empty()) {
        std::vector<std::string> pairs;
        for (const auto& [a, b] : replaced) pairs.push_back(a + " → " + b);
        parts.push_back("Default groups replaced by the chains of " + source + " (" + join(pairs) + "; " +
                        std::to_string(references.size()) + " reference(s) updated)");
    }
    if (!added.empty()) {
        parts.push_back("Added " + std::to_string(added.size()) + " group(s) from the chains of " + source + ": " +
                        join(added));
    }
    if (!pending.empty()) {
        parts.push_back(source + ": " + std::to_string(pending.size()) + " chain(s) without a group (" + join(pending) +
                        "); 'Groups from chains' adds them");
    }
    if (parts.empty()) {
        parts.push_back(!covered.empty() ? "Every chain of " + source + " already has a group"
                                         : "No polymer chains found in " + source);
    }
    std::string out;
    for (std::size_t k = 0; k < parts.size(); ++k) out += (k ? ". " : "") + parts[k];
    return out + ".";
}

std::vector<bool> polymerMask(const core::Topology& top) {
    const std::size_t n = top.residueCount();
    std::vector<bool> polymer(n, true);
    std::shared_ptr<const core::Topology> shared(&top, [](const core::Topology*) {});
    for (std::size_t k = 0; k < n; ++k) {
        const auto atoms = top.residueAtoms(static_cast<std::int32_t>(k));
        const std::string name = upper(strip(top.resname(static_cast<std::int32_t>(k))));
        bool water = kWaterResnames.count(name) > 0;
        const bool ion = atoms.size() == 1 || kIonResnames.count(name) > 0;
        if (atoms.size() <= 5 && !water && !ion) {  // water models by composition
            const core::AtomGroup group(shared, std::vector<std::int32_t>(atoms.begin(), atoms.end()));
            const std::vector<std::string> elements = core::atomElements(group);
            std::vector<std::string> heavy;
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const bool light = top.masses && (*top.masses)[static_cast<std::size_t>(atoms[i])] < 0.5;  // virtual sites
                if (elements[i] != "H" && !light) heavy.push_back(elements[i]);
            }
            water = heavy == std::vector<std::string>{"O"};
        }
        polymer[k] = !(water || ion);
    }
    return polymer;
}

std::vector<PolymerChain> polymerChains(const std::shared_ptr<const core::Topology>& topology,
                                        const std::vector<float>& positions) {
    const core::Topology& top = *topology;
    const std::vector<bool> mask = polymerMask(top);
    std::vector<std::int32_t> residues;
    for (std::size_t k = 0; k < mask.size(); ++k) {
        if (mask[k]) residues.push_back(static_cast<std::int32_t>(k));
    }
    if (residues.empty()) return {};
    std::vector<std::int32_t> firstAtoms;
    for (const std::int32_t r : residues) firstAtoms.push_back(top.residueAtoms(r).front());

    std::vector<std::string> chainIds;
    if (top.chainIDs) {
        for (const std::int32_t a : firstAtoms) chainIds.push_back(strip(top.strings[(*top.chainIDs)[static_cast<std::size_t>(a)]]));
    }
    std::vector<std::string> segids;
    for (const std::int32_t r : residues) segids.push_back(strip(top.segid(top.residueSegment[static_cast<std::size_t>(r)])));

    std::vector<std::string> labels;
    std::string source;
    if (usable(chainIds)) {
        labels = chainIds;
        source = "chainID";
    } else if (usable(segids)) {
        labels = segids;
        source = "segid";
    } else {
        // bonded fragments, else consecutive residues within covalent distance, else one fragment
        std::vector<std::int64_t> frag;
        if (top.bonds && !top.bonds->empty()) {
            const std::vector<std::int32_t> all = core::fragmentIndices(top);
            std::vector<std::int32_t> order;
            for (const std::int32_t a : firstAtoms) {
                const std::int32_t f = all[static_cast<std::size_t>(a)];
                auto it = std::find(order.begin(), order.end(), f);
                if (it == order.end()) {
                    order.push_back(f);
                    it = order.end() - 1;
                }
                frag.push_back(it - order.begin());
            }
        } else if (positions.size() == 3 * top.atomCount()) {
            std::shared_ptr<const core::Topology> shared = topology;
            std::vector<std::vector<std::int32_t>> heavyAtoms;
            for (const std::int32_t r : residues) {
                const auto atoms = top.residueAtoms(r);
                const core::AtomGroup group(shared, std::vector<std::int32_t>(atoms.begin(), atoms.end()));
                const std::vector<std::string> elements = core::atomElements(group);
                std::vector<std::int32_t> heavy;
                for (std::size_t i = 0; i < atoms.size(); ++i) {
                    if (elements[i] != "H") heavy.push_back(atoms[i]);
                }
                heavyAtoms.push_back(heavy.empty() ? std::vector<std::int32_t>(atoms.begin(), atoms.end()) : heavy);
            }
            std::int64_t k = 0;
            frag.push_back(0);
            for (std::size_t r = 0; r + 1 < heavyAtoms.size(); ++r) {
                float best = std::numeric_limits<float>::infinity();
                for (const std::int32_t a : heavyAtoms[r]) {
                    for (const std::int32_t b : heavyAtoms[r + 1]) {
                        float d2 = 0.0f;
                        for (std::size_t c = 0; c < 3; ++c) {
                            const float d = positions[3 * static_cast<std::size_t>(a) + c] -
                                            positions[3 * static_cast<std::size_t>(b) + c];
                            d2 += d * d;
                        }
                        best = std::min(best, std::sqrt(d2));
                    }
                }
                k += best <= kLinkMax ? 0 : 1;
                frag.push_back(k);
            }
        } else {
            frag.assign(residues.size(), 0);
        }
        for (const std::int64_t f : frag) labels.push_back(std::to_string(f + 1));
        source = "fragment";
        const auto single = [](const std::vector<std::string>& lab) {
            return !lab.empty() && std::all_of(lab.begin(), lab.end(), [](const std::string& s) { return !s.empty(); }) &&
                   std::set<std::string>(lab.begin(), lab.end()).size() == 1 && lab.front() != "SYSTEM";
        };
        if (*std::max_element(frag.begin(), frag.end()) == 0) {  // one chain: named after its chain ID / segment
            if (single(chainIds)) {
                labels = chainIds;
                source = "chainID";
            } else if (single(segids)) {
                labels = segids;
                source = "segid";
            }
        }
    }
    std::vector<PolymerChain> out;
    for (std::size_t k = 0; k < residues.size(); ++k) {
        auto it = std::find_if(out.begin(), out.end(), [&](const PolymerChain& c) { return c.label == labels[k]; });
        if (it == out.end()) {
            out.push_back({labels[k], source, {}});
            it = out.end() - 1;
        }
        it->residues.push_back(residues[k]);
    }
    return out;
}

bool untouchedDefaults(const std::vector<core::ResidueGroup>& groups) {
    const core::Project defaults;
    if (groups.size() != defaults.groups.size()) return false;
    for (std::size_t k = 0; k < groups.size(); ++k) {
        const core::ResidueGroup& a = groups[k];
        const core::ResidueGroup& b = defaults.groups[k];
        if (a.name != b.name || a.residues != b.residues || a.selection != b.selection || a.color != b.color ||
            !a.unknownKeys.empty()) {
            return false;
        }
    }
    return true;
}

std::vector<ChainGroup> chainGroups(const std::shared_ptr<const core::Topology>& topology, const core::Frame* current,
                                    const std::vector<core::ResidueGroup>& existing) {
    std::set<std::string> taken;
    for (const auto& g : existing) taken.insert(g.name);
    Colours colours(existing);
    std::map<std::vector<std::int32_t>, std::string> known;
    for (const auto& g : existing) {
        const auto ix = atomsOf(topology, current, g, existing);
        if (ix && !ix->empty()) known.emplace(*ix, g.name);  // the first group wins
    }
    static const std::vector<float> kNone;
    std::vector<ChainGroup> out;
    for (const PolymerChain& chain : polymerChains(topology, current ? current->positions : kNone)) {
        const std::string name = freeName("Chain_" + nameSafe(chain.label), taken);
        if (!core::groupNameProblem(name).empty()) continue;
        const std::vector<std::int32_t> want = atomsOfResidues(*topology, chain.residues);
        std::optional<core::ResidueGroup> found;
        for (const auto& [residues, selection] : specCandidates(*topology, chain)) {
            core::ResidueGroup group{name, residues, selection, ""};
            const auto got = atomsOf(topology, current, group, {group});
            if (got && *got == want) {
                found = group;
                break;
            }
        }
        if (!found) continue;  // an index selection always matches
        taken.insert(name);
        found->color = colours.next();
        const auto cover = known.find(want);
        out.push_back({*found, chain.label, chain.source, chain.residues.size(),
                       cover == known.end() ? std::string() : cover->second});
    }
    return out;
}

std::vector<ChainGroup> topologyChainGroups(const std::filesystem::path& path,
                                            const std::vector<core::ResidueGroup>& existing) {
    io::TopologyFile file = io::readTopologyFile(path, io::isParm7(core::pathToUtf8(path)) ? "PARM7" : "");
    const core::Frame* current = nullptr;
    if (file.coordinates && file.coordinates->frameCount() > 0) {
        current = &file.coordinates->seek(0);
    }
    return chainGroups(file.topology, current, existing);
}

ChainMerge mergeChainGroups(core::Project& project, const std::vector<ChainGroup>& chains, bool automatic) {
    ChainMerge out;
    if (chains.empty()) return out;
    const bool defaults = automatic && untouchedDefaults(project.groups);
    std::set<std::string> coveredBy;
    for (const auto& c : chains) coveredBy.insert(c.coveredBy);
    const bool fits = std::all_of(project.groups.begin(), project.groups.end(),
                                  [&](const core::ResidueGroup& g) { return coveredBy.count(g.name) > 0; });
    std::vector<const ChainGroup*> rest;
    if (defaults && !fits && chains.size() >= project.groups.size()) {
        const std::size_t n = project.groups.size();
        for (std::size_t k = 0; k < n; ++k) {
            const std::string old = project.groups[k].name;
            std::set<std::string> others;
            for (std::size_t j = 0; j < n; ++j) {
                if (j != k) others.insert(project.groups[j].name);
            }
            const std::string name = freeName(chains[k].group.name, others);
            for (std::string& ref : core::renameGroup(project, old, name)) {
                if (std::find(out.references.begin(), out.references.end(), ref) == out.references.end()) {
                    out.references.push_back(std::move(ref));
                }
            }
            project.groups[k].residues = chains[k].group.residues;  // keeps the default's colour
            project.groups[k].selection = chains[k].group.selection;
            out.replaced.emplace_back(old, name);
        }
        for (std::size_t k = n; k < chains.size(); ++k) rest.push_back(&chains[k]);
    } else if (project.groups.empty()) {
        for (const auto& c : chains) rest.push_back(&c);
    } else {
        for (const auto& c : chains) {
            if (!c.coveredBy.empty()) {
                out.covered.push_back(c.group.name);
            } else {
                rest.push_back(&c);
            }
        }
        if (automatic) {
            for (const ChainGroup* c : rest) out.pending.push_back(c->group.name);
            return out;
        }
    }
    for (const ChainGroup* c : rest) {
        const std::vector<std::string> names = core::groupNames(project);
        const std::string name = freeName(c->group.name, std::set<std::string>(names.begin(), names.end()));
        project.groups.push_back({name, c->group.residues, c->group.selection, c->group.color});
        out.added.push_back(name);
    }
    return out;
}

}  // namespace mdbridge::app
