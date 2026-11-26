// What a project cites: the settings' references and the run-time ones of its cached results.
#include "citations.hpp"

#include <algorithm>
#include <exception>
#include <optional>
#include <string_view>

#include "mdbridge/analysis/cache.hpp"

namespace mdbridge::app {

using core::Table;
using core::Value;

namespace {

bool contains(const std::vector<std::string>& items, const std::string& item) {
    return std::find(items.begin(), items.end(), item) != items.end();
}

std::string text(std::string_view view) { return std::string(view); }

core::Array textArray(const std::vector<std::string>& items) {
    core::Array out;
    for (const std::string& s : items) {
        out.emplace_back(s);
    }
    return out;
}

}  // namespace

ProjectCitations projectCitations(const core::Project& project) {
    namespace cit = core::citations;
    ProjectCitations out;
    const std::vector<std::string> resources = cit::projectResources(project, {});
    out.groups = cit::grouped(resources);
    out.settingsKeys = cit::keysFor(resources);
    try {
        const analysis::CacheStore cache(core::cacheDir(project));
        const std::optional<Value> manifest = cache.readManifest();
        const Value* entries = manifest && manifest->isTable() ? manifest->find("entries") : nullptr;
        if (entries && entries->isArray()) {
            for (const Value& entry : *entries->get<core::Array>()) {
                const Value* name = entry.isTable() ? entry.find("name") : nullptr;
                const Value* key = entry.isTable() ? entry.find("key") : nullptr;
                const std::optional<analysis::CachedResult> cached =
                    name && key && name->isString() && key->isString()
                        ? cache.loadResult(*name->get<std::string>(), *key->get<std::string>())
                        : std::nullopt;
                if (!cached) {  // an incomplete cache records nothing
                    out.cacheKeys.clear();
                    break;
                }
                for (const auto& [field, value] : cached->result.meta) {
                    if (field != "citation_keys" || !value.isArray()) {
                        continue;
                    }
                    for (const Value& item : *value.get<core::Array>()) {
                        if (const std::string* k = item.get<std::string>(); k && !contains(out.cacheKeys, *k)) {
                            out.cacheKeys.push_back(*k);
                        }
                    }
                }
            }
        }
    } catch (const std::exception&) {  // an unreadable manifest: no run-time keys
        out.cacheKeys.clear();
    }
    out.keys = cit::withRunKeys(out.settingsKeys, out.cacheKeys);
    for (const std::string& k : out.keys) {
        if (!contains(out.settingsKeys, k)) {
            out.runOnly.push_back(k);
        }
    }
    return out;
}

std::string citationText(const ProjectCitations& c) {
    if (c.keys.empty()) {
        return "This project uses no external resource with a reference yet.\n";
    }
    std::string out = "Resources this project uses (" + std::to_string(c.keys.size()) +
                      " references; each generated script lists only its own):\n\n";
    for (const auto& [res, refs] : c.groups) {
        if (refs.empty()) {
            continue;
        }
        out += text(res.name) + (res.note.empty() ? "" : " - " + text(res.note)) + "\n";
        for (const core::citations::Reference& r : refs) {
            out += "    " + text(r.text) + "\n";
        }
        out += "\n";
    }
    if (!c.runOnly.empty()) {
        out += "Recorded by the cached results of the last run:\n";
        for (const core::citations::Reference& r : core::citations::references(c.runOnly)) {
            out += "    " + text(r.text) + "\n";
        }
        out += "\n";
    }
    return out;
}

Value citationJson(const ProjectCitations& c) {
    core::Array resources;
    for (const auto& [res, refs] : c.groups) {
        std::vector<std::string> keys;
        for (const std::string_view k : res.keys) {
            keys.push_back(text(k));
        }
        resources.emplace_back(Table{{"id", text(res.id)},
                                     {"name", text(res.name)},
                                     {"note", text(res.note)},
                                     {"keys", textArray(keys)}});
    }
    return Table{{"keys", textArray(c.keys)},
                 {"settings_keys", textArray(c.settingsKeys)},
                 {"run_keys", textArray(c.cacheKeys)},
                 {"resources", std::move(resources)}};
}

}  // namespace mdbridge::app
