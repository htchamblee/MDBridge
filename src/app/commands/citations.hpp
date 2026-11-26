#pragma once
// What a project cites (owner decision D14, AP-16; docs/DIFFERENCES.md CMD-06): the references of
// the resources its settings use (core::citations::projectResources, keysFor) followed by the
// run-time references its cached results record (meta "citation_keys" of the entries of the cache
// manifest's latest run: the energy analysis' force fields, the covalent radii of a bond guess;
// core::citations::withRunKeys, CIT-06). Used by the `citations` command and the interface's About
// and Export citations (BibTeX) actions. No Qt.

#include <string>
#include <utility>
#include <vector>

#include "mdbridge/core/citations.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::app {

struct ProjectCitations {
    // the settings' resources with their references, in first-use order (core::citations::grouped)
    std::vector<std::pair<core::citations::Resource, std::vector<core::citations::Reference>>> groups;
    std::vector<std::string> settingsKeys;  // keysFor(the resources)
    std::vector<std::string> runKeys;       // keys the cached results add (not among settingsKeys)
    std::vector<std::string> runInstances;  // the cached results that record them (manifest order)
    std::vector<std::string> keys;          // withRunKeys(settingsKeys, all run keys): what is cited
};

// The project's citations. The cache is read from core::cacheDir(project) (project paths resolved);
// a missing or unreadable manifest or entry adds nothing.
ProjectCitations projectCitations(const core::Project& project);

// The text of the About dialog and of `citations` (MDATK's citation_text() for the settings, then
// the run-time part):
//   Resources this project uses (N references; each generated script lists only its own):
//   <blank>
//   <resource name>[ - <note>]
//       <reference text>
//   <blank>
//   Used at run time by the cached results (<instance>, ...):
//       <reference text>
//   <blank>
// N counts every key; without any: "This project uses no external resource with a reference yet."
// Lines end with "\n", the last one too.
std::string citationText(const ProjectCitations& citations);

// The --json form: {"keys": [...], "resources": [{"id", "name", "note", "keys"}], "run_keys":
// [...], "run_instances": [...], "references": [{"key", "text", "doi", "url", "bibtex"}]} with the
// references in `keys` order.
core::Value citationJson(const ProjectCitations& citations);

}  // namespace mdbridge::app
