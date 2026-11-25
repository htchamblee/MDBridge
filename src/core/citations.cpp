#include "mdbridge/core/citations.hpp"

#include <algorithm>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core::citations {

namespace {

#include "citation_data.inc"

void addUnique(std::vector<std::string>& out, std::string_view item) {
    if (!item.empty() && std::find(out.begin(), out.end(), item) == out.end()) {
        out.emplace_back(item);
    }
}

std::vector<std::string> wrap(std::string_view text, int width) {
    if (width <= 0) {
        throw ValueError("invalid width " + std::to_string(width) + " (must be > 0)");
    }
    return pythonTextwrap(text, width);
}

}  // namespace

std::span<const Reference> referenceTable() { return kReferences; }

std::span<const Resource> resourceTable() { return kResources; }

const Reference* findReference(std::string_view key) {
    for (const Reference& r : kReferences) {
        if (r.key == key) {
            return &r;
        }
    }
    return nullptr;
}

const Resource* findResource(std::string_view id) {
    for (const Resource& r : kResources) {
        if (r.id == id) {
            return &r;
        }
    }
    return nullptr;
}

std::vector<std::string> unique(std::span<const std::string> items) {
    std::vector<std::string> out;
    for (const std::string& item : items) {
        addUnique(out, item);
    }
    return out;
}

std::vector<std::string> keysFor(std::span<const std::string> resourceIds) {
    std::vector<std::string> out;
    for (const std::string& id : resourceIds) {
        if (const Resource* res = findResource(id)) {
            for (const std::string_view key : res->keys) {
                if (findReference(key)) {
                    addUnique(out, key);
                }
            }
        }
    }
    return out;
}

std::vector<Reference> references(std::span<const std::string> keys) {
    std::vector<Reference> out;
    for (const std::string& key : unique(keys)) {
        if (const Reference* r = findReference(key)) {
            out.push_back(*r);
        }
    }
    return out;
}

std::string bibtex(std::span<const std::string> keys) {
    std::string out;
    for (const Reference& r : references(keys)) {
        if (!out.empty()) {
            out += "\n\n";
        }
        out += pythonStrip(r.bibtex);
    }
    return out + "\n";
}

std::vector<std::string> textLines(std::span<const std::string> keys) {
    std::vector<std::string> out;
    for (const Reference& r : references(keys)) {
        out.emplace_back(r.text);
    }
    return out;
}

std::vector<std::pair<Resource, std::vector<Reference>>> grouped(
    std::span<const std::string> resourceIds) {
    std::vector<std::pair<Resource, std::vector<Reference>>> out;
    for (const std::string& id : unique(resourceIds)) {
        if (const Resource* res = findResource(id)) {
            std::vector<std::string> keys(res->keys.begin(), res->keys.end());
            out.emplace_back(*res, references(keys));
        }
    }
    return out;
}

std::vector<std::string> commentLines(std::span<const std::string> resourceIds,
                                      std::string_view prefix, int width,
                                      std::optional<std::string_view> heading) {
    const std::vector<Reference> refs = references(keysFor(resourceIds));
    if (refs.empty()) {
        return {};
    }
    const int room = width - static_cast<int>(characterCount(prefix));
    std::vector<std::string> lines;
    if (heading && !heading->empty()) {
        std::string names;
        for (const std::string& id : unique(resourceIds)) {
            if (const Resource* res = findResource(id)) {
                names += (names.empty() ? "" : ", ") + std::string(res->name);
            }
        }
        const std::string text = std::string(*heading) + " (" + names + "):";
        for (const std::string& w : wrap(text, room)) {
            lines.push_back(std::string(prefix) + w);
        }
    }
    for (const Reference& ref : refs) {
        std::vector<std::string> wrapped = wrap(ref.text, room - 4);
        if (wrapped.empty()) {
            wrapped.emplace_back(ref.text);
        }
        lines.push_back(std::string(prefix) + "  - " + wrapped[0]);
        for (std::size_t i = 1; i < wrapped.size(); ++i) {
            lines.push_back(std::string(prefix) + "    " + wrapped[i]);
        }
    }
    return lines;
}

std::vector<std::string> withRunKeys(std::span<const std::string> settingsKeys, std::span<const std::string> runKeys) {
    std::vector<std::string> out;
    for (const std::span<const std::string> part : {settingsKeys, runKeys}) {
        for (const std::string& key : part) {
            if (findReference(key) != nullptr && std::find(out.begin(), out.end(), key) == out.end()) {
                out.push_back(key);
            }
        }
    }
    return out;
}

std::vector<std::string> pipelineResources(const Project& project) {
    // CIT-01: MDATK starts with mdanalysis, numpy, scipy and matplotlib, the software it ran on;
    // MDBridge does not use them and cites the originators of the methods it uses instead.
    // CIT-02: never barnaba (MDBridge runs its builtin detector); RN-16: never vtk (MDBridge's
    // "vtk" is its Vulkan renderer).
    std::vector<std::string> out;
    if (project.cleaning.enabled) {
        addUnique(out, "cpptraj");
    }
    const RMSDParams& rmsd = project.analyses.rmsd;
    if (rmsd.enabled && !rmsd.targets.empty() && rmsd.superposition) {
        const auto fitted = [](const RMSDTarget& t) { return !pythonStrip(t.fit_selection).empty(); };
        if (!std::all_of(rmsd.targets.begin(), rmsd.targets.end(), fitted)) {
            addUnique(out, "qcp");  // superposedRmsd (rmsd_qcp.cpp)
        }
        if (std::any_of(rmsd.targets.begin(), rmsd.targets.end(), fitted)) {
            addUnique(out, "kabsch");  // the fit selection's superposition (spatial.cpp kabsch)
        }
    }
    if (project.output.video.enabled) {
        addUnique(out, "ffmpeg");
    }
    return out;
}

std::vector<std::string> projectResources(const Project& project,
                                          std::span<const std::vector<std::string>> stepResources) {
    std::vector<std::string> out;
    for (const std::vector<std::string>& step : stepResources) {
        for (const std::string& id : step) {
            addUnique(out, id);
        }
    }
    for (const std::string& id : pipelineResources(project)) {
        addUnique(out, id);
    }
    return out;
}

}  // namespace mdbridge::core::citations
