#pragma once
// Problems that say where they are: the presentation of MDATK 0.0.10 (mdatk/problems.py of MDATK
// 0.1.1) over MDBridge's checks.
//
// What is a problem stays what validate() reports (MDATK v0.0.4.1's blocking rules with
// MDBridge's fixes, owner decision D7 of milestone 5): projectProblems() gives one located
// Problem per line of validate(), in the same order, followed by the warnings. Each carries the
// project path of its setting (the field paths of project_document.hpp), the GUI page that shows
// it and a readable location ("Analyses › Hydrogen bonds › Group A"); its message is MDATK 0.1.1's
// text for the same check where 0.1.1 has one, else validate()'s.

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::core {

struct Project;

// The GUI side-bar pages, in MDATK's order: id and title.
struct ProblemPage {
    std::string_view id;
    std::string_view title;
};
inline constexpr ProblemPage kProblemPages[] = {
    {"input", "Input & cleaning"}, {"groups", "Residue groups"}, {"analyses", "Analyses"},
    {"structure", "Structure view"}, {"layout", "Layout & panels"}, {"outputs", "Outputs"},
};

struct Problem {
    // Field path of the setting ("analyses.hbonds.group_a", "groups.1.residues", "layout.panels";
    // "" when the problem names no setting).
    std::string field;
    // What is wrong, without the location (MDATK 0.1.1's words where it has the same check).
    std::string message;
    // The line validate() gives for it (the command line prints these; "" for a warning).
    std::string check;
    // "error" (validate() reports it: a run is refused) or "warning".
    std::string level = "error";
    // describeField(field, project): "Residue groups › Group 2 (Target) › Residues".
    std::string location;

    // The page that shows the field: problemPage(field).
    std::string_view page() const;
    // "<location>: <message>" (the message alone without a location), as MDATK 0.1.1 shows it.
    std::string text() const;

    bool operator==(const Problem&) const = default;
};

// tab_of: the page id of a field path ("" for a path no page shows, e.g. "workflow.nodes").
std::string_view problemPage(std::string_view field);

// The title of a page id ("" for an unknown id).
std::string_view problemPageTitle(std::string_view page);

// describe(field, project) of MDATK 0.1.1: the page title (or the first component, capitalised),
// then a label per known level of the path; list items are numbered from 1 and, with a project,
// named after the item ("Group 2 (Target)", panels after their type: "Panel 3 (Contact map)"); the
// last component gets a leaf label or its own name, capitalised with '_' as spaces. "" for "".
// "layout.panels.2.groups" -> "Layout & panels › Panel 3 (Contact map) › Groups".
std::string describeField(std::string_view field, const Project* project = nullptr);

// Every problem of the project: the errors (one per validate() line, same order; Problem::check
// is that line) and then the warnings (layoutWarnings() at "layout.panels"). Files are looked up
// relative to the current directory, as validate() does: pass a project with resolved paths
// (ProjectDocument::problems does).
std::vector<Problem> projectProblems(const Project& project);

// A problem of `field` with its location described from `project` (may be null).
Problem locatedProblem(std::string field, std::string message, const Project* project,
                       std::string level = "error");

}  // namespace mdbridge::core
