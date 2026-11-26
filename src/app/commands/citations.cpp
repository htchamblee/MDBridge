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

    std::vector<std::string> recorded;  // every run-time key, in manifest order
    try {
        const analysis::CacheStore cache(core::cacheDir(project));
        const std::optional<Value> manifest = cache.readManifest();
        const Value* entries = manifest && manifest->isTable() ? manifest->find("entries") : nullptr;
        if (entries && entries->isArray()) {
            for (const Value& entry : *entries->get<core::Array>()) {
                const Value* name = entry.isTable() ? entry.find("name") : nullptr;
                const Value* key = entry.isTable() ? entry.find("key") : nullptr;
                if (!name || !key || !name->isString() || !key->isString()) {
                    continue;
                }
                const std::optional<analysis::CachedResult> cached =
                    cache.loadResult(*name->get<std::string>(), *key->get<std::string>());
                if (!cached) {
                    continue;
                }
                bool adds = false;
                for (const auto& [field, value] : cached->result.meta) {
                    if (field != "citation_keys" || !value.isArray()) {
                        continue;
                    }
                    for (const Value& item : *value.get<core::Array>()) {
                        const std::string* k = item.get<std::string>();
                        if (k == nullptr || cit::findReference(*k) == nullptr) {
                            continue;
                        }
                        if (!contains(recorded, *k)) {
                            recorded.push_back(*k);
                        }
                        if (!contains(out.settingsKeys, *k)) {
                            adds = true;
                            if (!contains(out.runKeys, *k)) {
                                out.runKeys.push_back(*k);
                            }
                        }
                    }
                }
                if (adds && !contains(out.runInstances, *name->get<std::string>())) {
                    out.runInstances.push_back(*name->get<std::string>());
                }
            }
        }
    } catch (const std::exception&) {  // an unreadable manifest: no run-time keys
    }
    out.keys = cit::withRunKeys(out.settingsKeys, recorded);
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
    if (!c.runKeys.empty()) {
        std::string instances;
        for (std::size_t k = 0; k < c.runInstances.size(); ++k) {
            instances += (k ? ", " : "") + c.runInstances[k];
        }
        out += "Used at run time by the cached results (" + instances + "):\n";
        for (const core::citations::Reference& r : core::citations::references(c.runKeys)) {
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
        for (const core::citations::Reference& r : refs) {
            keys.push_back(text(r.key));
        }
        resources.emplace_back(Table{{"id", text(res.id)},
                                     {"name", text(res.name)},
                                     {"note", text(res.note)},
                                     {"keys", textArray(keys)}});
    }
    core::Array references;
    for (const core::citations::Reference& r : core::citations::references(c.keys)) {
        references.emplace_back(Table{{"key", text(r.key)},
                                      {"text", text(r.text)},
                                      {"doi", text(r.doi)},
                                      {"url", text(r.url)},
                                      {"bibtex", text(r.bibtex)}});
    }
    return Table{{"keys", textArray(c.keys)},
                 {"resources", std::move(resources)},
                 {"run_keys", textArray(c.runKeys)},
                 {"run_instances", textArray(c.runInstances)},
                 {"references", std::move(references)}};
}

}  // namespace mdbridge::app
