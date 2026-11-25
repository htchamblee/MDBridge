// Instance ids, data keys, labels and the residue groups of cache keys (mdatk/analysis/base.py,
// mdatk/analysis/__init__.py expected_keys).

#include <map>
#include <set>
#include <stdexcept>

#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::analysis {

namespace {

std::string repr(std::string_view text) {
    return core::pythonStringRepr(std::string(text));
}

// [A-Za-z0-9_.+\-'] of base.py's _GROUP_TOKEN
bool isGroupTokenCharacter(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.' ||
           c == '+' || c == '-' || c == '\'';
}

// _GROUP_TOKEN.findall(text): the names after every "group:"
std::vector<std::string> groupTokens(std::string_view text) {
    static constexpr std::string_view kPrefix = "group:";
    std::vector<std::string> out;
    std::size_t pos = text.find(kPrefix);
    while (pos != std::string_view::npos) {
        std::size_t end = pos + kPrefix.size();
        while (end < text.size() && isGroupTokenCharacter(text[end])) {
            ++end;
        }
        if (end > pos + kPrefix.size()) {
            out.emplace_back(text.substr(pos + kPrefix.size(), end - pos - kPrefix.size()));
        }
        pos = text.find(kPrefix, end > pos + kPrefix.size() ? end : pos + 1);
    }
    return out;
}

std::string joined(const std::vector<std::string>& items) {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        out += (i ? ", " : "") + items[i];
    }
    return out;
}

}  // namespace

Instance parseInstance(std::string_view instance) {
    Instance out;
    const std::size_t at = instance.find(kInstanceSeparator);
    out.base = std::string(instance.substr(0, at));
    if (at == std::string_view::npos) {
        return out;
    }
    std::string_view rest = instance.substr(at + 1);
    while (!rest.empty()) {
        const std::size_t plus = rest.find(kGroupSeparator);
        const std::string_view group = rest.substr(0, plus);
        if (!group.empty()) {
            out.groups.emplace_back(group);
        }
        if (plus == std::string_view::npos) {
            break;
        }
        rest = rest.substr(plus + 1);
    }
    return out;
}

std::string instanceId(std::string_view base, const std::vector<std::string>& groups) {
    std::vector<std::string> names;
    for (const std::string& g : groups) {
        if (!g.empty()) {
            names.push_back(g);
        }
    }
    std::string out(base);
    if (!names.empty()) {
        out += kInstanceSeparator;
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (i) {
                out += kGroupSeparator;
            }
            out += names[i];
        }
    }
    return out;
}

std::pair<std::string, std::string> splitKey(std::string_view key) {
    const std::size_t slash = key.find('/');
    if (slash == std::string_view::npos) {
        return {std::string(key), {}};
    }
    return {std::string(key.substr(0, slash)), std::string(key.substr(slash + 1))};
}

const AnalysisDescriptor& analysisFor(std::string_view instance) {
    const std::string base = parseInstance(instance).base;
    if (const AnalysisDescriptor* d = findAnalysis(base)) {
        return *d;
    }
    std::vector<std::string> known;
    for (const AnalysisDescriptor* d : registeredAnalyses()) {
        known.push_back(d->name());
    }
    throw std::out_of_range("unknown analysis " + repr(base) + " (known: " + joined(known) + ")");
}

core::Project instanceProject(const core::Project& project, std::string_view instance) {
    const AnalysisDescriptor& d = analysisFor(instance);
    const std::vector<std::string> groups = parseInstance(instance).groups;
    if (groups.empty()) {
        return project;
    }
    if (d.groupSlots() == 0) {
        throw core::ValueError(d.name() + " does not take residue groups (" + repr(instance) + ")");
    }
    if (groups.size() != static_cast<std::size_t>(d.groupSlots())) {
        throw core::ValueError(repr(instance) + ": " + d.name() + " takes " + std::to_string(d.groupSlots()) +
                               " group(s)");
    }
    const std::vector<std::string> names = core::groupNames(project);
    const std::set<std::string> known(names.begin(), names.end());
    std::vector<std::string> missing;
    for (const std::string& g : groups) {
        if (!known.contains(g)) {
            missing.push_back(g);
        }
    }
    if (!missing.empty()) {
        throw core::ValueError(repr(instance) + ": unknown residue group(s) " + joined(missing));
    }
    core::Project copy = project;
    d.withGroups(copy, groups);
    return copy;
}

std::string keyLabel(std::string_view key, const core::Project& project) {
    const auto [instance, out] = splitKey(key);
    try {
        const AnalysisDescriptor& d = analysisFor(instance);
        return d.outputLabel(instanceProject(project, instance), out);
    } catch (const std::out_of_range&) {
        return std::string(key);
    } catch (const core::ValueError&) {
        return std::string(key);
    }
}

ExpectedOutputs expectedKeys(const core::Project& project, bool includeDisabled) {
    ExpectedOutputs out;
    for (const AnalysisDescriptor* d : registeredAnalyses()) {
        if (!includeDisabled && !d->internal() && !d->enabled(project)) {
            continue;
        }
        const ExpectedOutputs outputs = d->expectedOutputs(project);
        for (const std::string& n : outputs.series) {
            out.series.push_back(d->name() + "/" + n);
        }
        for (const std::string& n : outputs.matrices) {
            out.matrices.push_back(d->name() + "/" + n);
        }
        for (const std::string& n : outputs.pairmaps) {
            out.pairmaps.push_back(d->name() + "/" + n);
        }
    }
    return out;
}

core::Array referencedGroups(const std::vector<std::string>& texts, const core::Project& project) {
    // a repeated name: its first group, the one selections use (MDATK: the last)
    std::map<std::string, const core::ResidueGroup*> byName;
    for (const core::ResidueGroup& g : project.groups) {
        byName.emplace(g.name, &g);
    }
    std::map<std::string, core::Table> seen;  // sorted by name, as sorted() of Python strings
    std::vector<std::string> stack(texts.begin(), texts.end());
    while (!stack.empty()) {
        const std::string text = core::pythonStrip(stack.back());
        stack.pop_back();
        if (text.empty()) {
            continue;
        }
        std::vector<std::string> names = groupTokens(text);
        if (byName.contains(text)) {
            names.push_back(text);
        }
        for (const std::string& name : names) {
            const auto group = byName.find(name);
            if (seen.contains(name) || group == byName.end()) {
                continue;
            }
            const core::ResidueGroup& g = *group->second;
            seen[name] = core::Table{{"name", g.name}, {"residues", g.residues}, {"selection", g.selection}};
            stack.push_back(g.selection);
        }
    }
    core::Array out;
    for (auto& [name, definition] : seen) {
        out.emplace_back(std::move(definition));
    }
    return out;
}

}  // namespace mdbridge::analysis
