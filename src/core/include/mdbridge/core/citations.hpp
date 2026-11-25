#pragma once
// Citations of the programs, force fields, parameter sets and methods a project uses
// (mdatk/citations.py).
//
// Only what a project actually uses is cited: each planner step reports the resources its
// settings turn on, the analysis and video pipeline reports its own (pipelineResources), and the
// resource table maps every resource to the references that describe it. A generated script cites
// its own step; a project (About dialog, BibTeX export, parameter log) cites every step plus the
// pipeline. The tables are MDATK's (src/core/citation_data.inc, tools/gen_citation_data.py) with
// MDBridge's changes (tests/parity/improvements_citations.py; docs/DIFFERENCES.md "Citations"):
//
// * Methods are cited by their originators, not by one implementation of them (owner decision).
//   MDATK cites MDAnalysis, NumPy, SciPy and Matplotlib for every project because it ran on them;
//   MDBridge does not use those packages and does not cite them (CIT-01). Where a method with a
//   published originator decides a reported result, that work is cited when the project uses it:
//   the QCP superposition (Theobald 2005; Liu et al. 2010) and Kabsch's (1976, 1978) for RMSD.
// * Standard algorithms and common knowledge (distance searches, minimum images, geometric
//   H-bond / contact / stacking criteria, Gaussian filters, linear interpolation, tick choice) are
//   not cited, nor are colour maps (styling) or the random numbers of the example data.
// * Data a result depends on is credited where it is used: the covalent radii of the bond guess
//   (Cordero et al. 2008) in the energy and atoms analyses' guessed-bonds warnings (CIT-03,
//   CIT-05). The energy analysis' parameters without parm7 (Amber force-field residue templates and
//   atom types) are cited by the run that uses them, in its log and metadata (CIT-04): which force
//   fields a run uses depends on its residues, not on the settings, so they are no project
//   resource; their references are in the table.

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

// REFERENCE_DATA / RESOURCE_DATA in MDATK's order, with MDBridge's changes (removals; additions at
// the end).
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

// Resources of the analysis, rendering and video stages: cpptraj (cleaning.enabled); qcp and
// kabsch (analyses.rmsd enabled with targets and superposition: qcp when a target has a blank
// fit_selection, kabsch when one has a fit selection); ffmpeg (output.video.enabled). Unlike
// MDATK: never mdanalysis, numpy, scipy or matplotlib (CIT-01), barnaba (MDBridge never runs it,
// CIT-02) or vtk (MDBridge's "vtk" is its Vulkan renderer, RN-16).
std::vector<std::string> pipelineResources(const Project& project);

// project_resources(project): everything the project uses, its planner steps and the pipeline.
// MDATK resolves the project's plan itself (Plan(workflow, workflow_dir, groups); a plan that
// cannot be resolved counts as none) and asks each step's tool for its resources. The planner
// is milestone 7: its caller passes, in plan order, the resources (tool.resources(settings)) of
// every resolved step whose tool generates a script; an empty list when there is no plan.
std::vector<std::string> projectResources(const Project& project,
                                          std::span<const std::vector<std::string>> stepResources);

}  // namespace mdbridge::core::citations
