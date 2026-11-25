#pragma once
// Citations of the programs, force fields and parameter sets a project uses (mdatk/citations.py).
//
// Only what a project actually uses is cited: each planner step reports the resources its
// settings turn on, the analysis and video pipeline reports its own (pipelineResources), and the
// resource table maps every resource to the references that describe it. A generated script cites
// its own step; a project (About dialog, BibTeX export, parameter log) cites every step plus the
// pipeline. The tables are MDATK's (src/core/citation_data.inc, tools/gen_citation_data.py).

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/project.hpp"

namespace mdbridge::core::citations {

struct Reference {
    std::string_view key;       // BibTeX key
    std::string_view text;      // one-line citation: authors, journal year, volume, pages, doi
    std::string_view bibtex;
    std::string_view doi;
    std::string_view url;
    std::string_view verified;  // where the details were checked: pdf:<file>, repo:<url>, ...
};

struct Resource {
    std::string_view id;
    std::string_view name;
    std::span<const std::string_view> keys;  // references to cite (in this order)
    std::string_view note;  // e.g. "cite together with ..." or how to cite software without a paper
};

// REFERENCE_DATA / RESOURCE_DATA in MDATK's order.
std::span<const Reference> referenceTable();
std::span<const Resource> resourceTable();
// REFERENCES[key] / RESOURCES[id], or nullptr.
const Reference* findReference(std::string_view key);
const Resource* findResource(std::string_view id);

// The non-empty items in first-occurrence order.
std::vector<std::string> unique(std::span<const std::string> items);

// Reference keys of the resources in first-use order (unknown ids and keys are ignored).
// step_citations(step) of MDATK is keysFor(<the step's resources>).
std::vector<std::string> keysFor(std::span<const std::string> resourceIds);

// The references of the keys (unique, unknown keys dropped).
std::vector<Reference> references(std::span<const std::string> keys);

// The BibTeX entries of the keys, stripped, separated by blank lines, ending with "\n".
std::string bibtex(std::span<const std::string> keys);

// The one-line citations of the keys.
std::vector<std::string> textLines(std::span<const std::string> keys);

// (resource, its references) for display, in first-use order; unknown ids are dropped.
std::vector<std::pair<Resource, std::vector<Reference>>> grouped(
    std::span<const std::string> resourceIds);

// The citation block of a script header: "<prefix><heading> (<resource names>):" wrapped to
// width - len(prefix) (no heading when it is unset or empty), then one entry per reference,
// "<prefix>  - " + its text wrapped to width - len(prefix) - 4, continued with "<prefix>    ".
// Empty when the resources have no references. Throws ValueError("invalid width N (must be > 0)")
// as textwrap does when a wrap width is not positive.
std::vector<std::string> commentLines(std::span<const std::string> resourceIds,
                                      std::string_view prefix = "# ", int width = 100,
                                      std::optional<std::string_view> heading = "Please cite");

// Resources of MDATK's own analysis, rendering and video stages: mdanalysis, numpy, scipy,
// matplotlib; cpptraj (cleaning.enabled), ffmpeg (output.video.enabled), barnaba
// (analyses.basepairs.engine "barnaba" and a secondary panel whose "source" option is not
// "dotbracket"). Not vtk: MDATK cites it for visual.renderer "vtk", which in MDBridge means its
// own GPU renderer (Vulkan; RN-16, docs/DIFFERENCES.md).
std::vector<std::string> pipelineResources(const Project& project);

// project_resources(project): everything the project uses, its planner steps and the pipeline.
// MDATK resolves the project's plan itself (Plan(workflow, workflow_dir, groups); a plan that
// cannot be resolved counts as none) and asks each step's tool for its resources. The planner
// is milestone 7: its caller passes, in plan order, the resources (tool.resources(settings)) of
// every resolved step whose tool generates a script; an empty list when there is no plan.
std::vector<std::string> projectResources(const Project& project,
                                          std::span<const std::vector<std::string>> stepResources);

}  // namespace mdbridge::core::citations
