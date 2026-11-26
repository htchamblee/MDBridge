// The GUI's edit model (core/project_document.hpp, problems.hpp; milestone 5 W0-B): round trips of
// every project file MDATK wrote, byte for byte, through open -> save; the lossless rules of AP-01 ..
// AP-06; get / set / text by field path; undo / redo; the group rename cascade and the located
// problems against MDATK 0.1.1 (data/edit-model, tests/parity/make_edit_model_cases.py); recovery
// copies.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/step_log.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const fs::path kProjects = kData / "projects";
const fs::path kGolden = kData / "project-file" / "golden";
const fs::path kWorkflow = kData / "project-file" / "workflow";
const fs::path kLatest = kData / "project-file" / "mdatk-0.1.1";
const fs::path kCases = kData / "edit-model" / "edit_model_cases.json";

std::string unixLines(std::string text) {
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '\r' && i + 1 < text.size() && text[i + 1] == '\n') {
            continue;
        }
        out += text[i];
    }
    return out;
}

std::string utf8(const fs::path& path) { return pathToUtf8(path); }

// "" when equal, else the first line that differs (for readable failures of long texts)
std::string firstDifference(const std::string& a, const std::string& b) {
    if (a == b) {
        return {};
    }
    std::size_t line = 1;
    std::size_t start = 0;
    for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i) {
        if (a[i] != b[i]) {
            const std::size_t endA = a.find('\n', i);
            const std::size_t endB = b.find('\n', i);
            return "line " + std::to_string(line) + ":\n  " + a.substr(start, endA - start) + "\n  " +
                   b.substr(start, endB - start);
        }
        if (a[i] == '\n') {
            ++line;
            start = i + 1;
        }
    }
    return "one text ends at line " + std::to_string(line) + " (sizes " + std::to_string(a.size()) + ", " +
           std::to_string(b.size()) + ")";
}

std::string plainText(const Project& p) { return toJson(toPlain(p), 1); }

// A clock that always says `now` (the `saved` stamp a file had), restored at the end.
class FixedClock {
public:
    explicit FixedClock(DateTime now) : m_previous(setClock(Clock{[now] { return now; }, {}})) {}
    ~FixedClock() { setClock(m_previous); }
    FixedClock(const FixedClock&) = delete;
    FixedClock& operator=(const FixedClock&) = delete;

private:
    Clock m_previous;
};

class ScopedCwd {
public:
    explicit ScopedCwd(const fs::path& path) : m_old(fs::current_path()) { fs::current_path(path); }
    ~ScopedCwd() { fs::current_path(m_old); }
    ScopedCwd(const ScopedCwd&) = delete;
    ScopedCwd& operator=(const ScopedCwd&) = delete;

private:
    fs::path m_old;
};

DateTime savedStamp(const fs::path& path) {
    const Value* saved = readProject(path).header.find("saved");
    return saved && saved->get<DateTime>() ? *saved->get<DateTime>() : DateTime{};
}

std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    if (const Array* array = value.get<Array>()) {
        for (const Value& item : *array) {
            out.push_back(item.stringOr());
        }
    }
    return out;
}

const Array& items(const Value& value, std::string_view key) { return *value.find(key)->get<Array>(); }

// MDATK 0.1.1's locations "groups[1].selection" as field paths "groups.1.selection"
std::vector<std::string> dotted(const std::vector<std::string>& paths) {
    static const std::regex index(R"(\[(\d+)\])");
    std::vector<std::string> out;
    for (const std::string& p : paths) {
        out.push_back(std::regex_replace(p, index, ".$1"));
    }
    return out;
}

// The project files MDATK wrote (folder, file): every one must come back byte for byte.
std::vector<std::pair<fs::path, std::string>> mdatkWrittenProjects() {
    std::vector<std::pair<fs::path, std::string>> out;
    for (const char* name : {"default.mdatk.toml", "full.mdatk.toml", "split.mdatk.toml"}) {
        out.emplace_back(kProjects, name);
    }
    for (const char* name : {"default.mdatk.toml", "rich.mdatk.toml", "rich.resaved.mdatk.toml",
                             "rich.with-log.mdatk.toml", "split-nolog.after.mdatk.toml", "split.after.mdatk.toml",
                             "split.before.mdatk.toml", "split.resaved.mdatk.toml", "truncation.mdatk.toml"}) {
        out.emplace_back(kGolden, name);
    }
    for (const auto& entry : fs::directory_iterator(kWorkflow)) {
        const std::string name = utf8(entry.path().filename());
        if (name.size() > 11 && name.substr(name.size() - 11) == ".mdatk.toml") {
            out.emplace_back(kWorkflow, name);
        }
    }
    for (const char* name : {"plan.mdatk.toml", "sidecar.mdatk.toml"}) {
        out.emplace_back(kLatest, name);
    }
    return out;
}

void copyFolder(const fs::path& from, const fs::path& to) {
    fs::create_directories(to);
    for (const auto& entry : fs::directory_iterator(from)) {
        if (entry.is_regular_file()) {
            fs::copy_file(entry.path(), to / entry.path().filename());
        }
    }
}

struct Events {
    std::vector<DocumentChange> seen;
    ProjectDocument::ListenerId id = 0;
    ProjectDocument* doc = nullptr;

    explicit Events(ProjectDocument& d) : doc(&d) {
        id = d.addListener([this](const DocumentChange& c) { seen.push_back(c); });
    }
    ~Events() { doc->removeListener(id); }
    Events(const Events&) = delete;
    Events& operator=(const Events&) = delete;
};

}  // namespace

// ---- round trips -----------------------------------------------------------------------------

// Every project file MDATK wrote (0.0.4.1 and 0.1.1) comes back byte for byte through open ->
// save without an edit, step logs and log files included; the paths stay as written and the
// resolved project is what loadProject reads.
TEST(ProjectDocument, EveryMdatkProjectRoundTripsByteForByte) {
    std::size_t count = 0;
    for (const auto& [folder, name] : mdatkWrittenProjects()) {
        TempDir dir;
        copyFolder(folder, dir.path());
        const fs::path path = dir / name;
        const std::string before = readFile(path);
        std::vector<std::string> logs;
        for (const auto& entry : fs::directory_iterator(dir.path())) {
            logs.push_back(utf8(entry.path().filename()));
        }
        ProjectDocument doc;
        doc.open(path);
        EXPECT_FALSE(doc.isModified()) << name;
        EXPECT_EQ(firstDifference(plainText(doc.resolvedProject()), plainText(loadProject(path))), "") << name;
        {
            const FixedClock clock(savedStamp(path));
            doc.save();
        }
        std::string expected = unixLines(before);
        if (name == "split.after.mdatk.toml" || name == "sidecar.mdatk.toml") {
            // MDATK's log move leaves a blank line at the end, which its own next save drops too
            // (golden/split.resaved.mdatk.toml)
            ASSERT_EQ(expected.substr(expected.size() - 2), "\n\n");
            expected.pop_back();
        }
        EXPECT_EQ(firstDifference(unixLines(readFile(path)), expected), "") << folder / name;
        for (const auto& entry : fs::directory_iterator(dir.path())) {  // side files untouched
            if (entry.path().filename() != fs::path(name) && entry.path().extension() == ".toml") {
                EXPECT_EQ(readFile(entry.path()), readFile(folder / entry.path().filename())) << entry.path();
            }
        }
        EXPECT_FALSE(doc.isModified()) << name;
        EXPECT_FALSE(doc.changedOnDisk()) << name;
        ++count;
    }
    EXPECT_GE(count, 18U);
}

// A file of a later MDATK keeps its comments (it says that MDATK wrote it, PRJ-09): an edit changes
// only its own line; an entry the file did not have gets MDBridge's comment.
TEST(ProjectDocument, LaterMdatkCommentsAreKept) {
    TempDir dir;
    copyFolder(kLatest, dir.path());
    const fs::path path = dir / "plan.mdatk.toml";
    const std::string before = unixLines(readFile(path));
    ProjectDocument doc;
    doc.open(path);
    EXPECT_TRUE(doc.comments());
    doc.setText("analyses.contacts.cutoff", "5");
    const FixedClock clock(savedStamp(path));
    doc.save();
    std::string expected = before;
    const std::string line = "\ncutoff = 4.5\n# Ignore hydrogens";
    ASSERT_NE(expected.find(line), std::string::npos);
    expected.replace(expected.find(line), line.size(), "\ncutoff = 5.0\n# Ignore hydrogens");
    EXPECT_EQ(firstDifference(unixLines(readFile(path)), expected), "");
    // a recovery copy keeps them too
    const fs::path copy = writeRecovery(doc, dir / "recovery");
    ProjectDocument restored;
    openRecovery(restored, findRecoveries(dir / "recovery").front());
    EXPECT_EQ(fs::path(findRecoveries(dir / "recovery").front().file), copy);
    restored.save();
    EXPECT_EQ(firstDifference(unixLines(readFile(path)), expected), "");
    // a 0.0.4.1 file gets MDBridge's comments
    ProjectDocument old;
    old.open(kProjects / "default.mdatk.toml");
    EXPECT_FALSE(old.comments());
}

// Every field set to the text it shows changes nothing: no undo step, not modified, and the file
// saved afterwards is the same (AP-01 .. AP-05: what a widget shows and gives back is lossless).
TEST(ProjectDocument, EveryFieldGivenItsOwnTextChangesNothing) {
    for (const auto& [folder, name] : mdatkWrittenProjects()) {
        ProjectDocument doc;
        doc.open(folder / name);
        const Project before = doc.project();
        std::size_t texts = 0;
        for (const std::string& path : fieldPaths(doc.project())) {
            const auto info = doc.info(path);
            ASSERT_TRUE(info) << name << " " << path;
            if (info->kind == FieldKind::NumberRows || info->kind == FieldKind::Records ||
                info->kind == FieldKind::Record) {
                EXPECT_THROW(doc.text(path), FieldError) << path;
                continue;
            }
            EXPECT_FALSE(doc.setText(path, doc.text(path))) << name << " " << path;
            // the value the text means is the value held (no rounding, no clamping)
            const Value parsed = parseFieldText(doc.project(), path, doc.text(path));
            EXPECT_FALSE(doc.set(path, parsed)) << name << " " << path << " = " << doc.text(path);
            ++texts;
        }
        EXPECT_GT(texts, 150U) << name;
        EXPECT_FALSE(doc.canUndo()) << name;
        EXPECT_FALSE(doc.isModified()) << name;
        EXPECT_EQ(plainText(doc.project()), plainText(before)) << name;
    }
}

// Files written by hand (aliases, CRLF, no marker) and parameter logs keep their values through a
// save; the second save is the first one's bytes.
TEST(ProjectDocument, HandMadeProjectsKeepTheirValues) {
    for (const char* name : {"handwritten.mdatk.toml", "datetimes.mdatk.toml", "missing_log.mdatk.toml",
                             "empty_params_log.toml", "MDATK_params_log.toml"}) {
        TempDir dir;
        fs::copy_file(kProjects / name, dir / name);
        ProjectDocument doc;
        doc.open(dir / name);
        const Value plain = toPlain(doc.project());
        const fs::path target = dir / "saved.mdatk.toml";
        if (doc.isParamsLog()) {
            EXPECT_THROW(doc.save(), ProjectError) << name;
            EXPECT_TRUE(doc.isParamsLog());
        }
        doc.saveAs(target);
        EXPECT_FALSE(doc.isParamsLog());
        EXPECT_EQ(doc.path(), pathFromUtf8(utf8(target))) << name;
        ProjectDocument again;
        again.open(target);
        EXPECT_EQ(pythonRepr(toPlain(again.project())), pythonRepr(plain)) << name;
        const std::string first = readFile(target);
        const FixedClock clock(savedStamp(target));
        again.save();
        EXPECT_EQ(readFile(target), first) << name;
    }
}

// AP-06 / D8: relative paths are written as they were, not made absolute.
TEST(ProjectDocument, PathsStayAsWritten) {
    TempDir dir;
    fs::copy_file(kProjects / "handwritten.mdatk.toml", dir / "p.mdatk.toml");
    ProjectDocument doc;
    doc.open(dir / "p.mdatk.toml");
    EXPECT_EQ(doc.project().inputs.topology, "system.parm7");
    EXPECT_EQ(doc.project().output.out_dir, "results");
    const Project resolved = doc.resolvedProject();
    EXPECT_EQ(resolved.inputs.topology, resolvePath("system.parm7", doc.folder()));
    EXPECT_EQ(resolved, loadProject(dir / "p.mdatk.toml"));
    doc.setText("inputs.trajectories", "md.nc, md2.nc");
    doc.save();
    const std::string text = readFile(dir / "p.mdatk.toml");
    EXPECT_NE(text.find("topology = \"system.parm7\""), std::string::npos);
    EXPECT_NE(text.find("trajectories = [\"md.nc\", \"md2.nc\"]"), std::string::npos);
    EXPECT_NE(text.find("out_dir = \"results\""), std::string::npos);
}

// The values MDATK's GUI rewrote on save (AP-01 .. AP-05) are saved back exactly.
TEST(ProjectDocument, ValuesMdatksWidgetsRewroteAreKept) {
    TempDir dir;
    Project p;
    p.analyses.contacts.cutoff = 25.0;      // AP-01: spin range 20
    p.analyses.hbonds.angle_cutoff = 60.0;  // 90
    p.analyses.energy.cutoff = 60.0;        // 50
    p.output.video.width = 7682;            // 7680
    p.visual.zoom = 0.01;                   // 0.05
    p.inputs.dt_ps = 0.0123456789;          // AP-02: %g gave 0.0123457
    p.analyses.energy.target_value = -1234567.5;
    p.visual.base_color = "tab:orange";     // AP-03: GDK made these black or #rrggbb
    p.visual.glow_color = "#ff000080";
    p.groups[0].color = "gold";
    p.output.charts.formats = {"jpg", "svg"};  // AP-04: check boxes gave ["svg"]
    p.visual.renderer = "vulkan";           // AP-05: combos replaced unknown values
    p.visual.color_mode = "rainbow";
    p.inputs.time_unit = "fs";
    saveProject(dir / "p.mdatk.toml", p);
    const std::string before = readFile(dir / "p.mdatk.toml");

    ProjectDocument doc;
    doc.open(dir / "p.mdatk.toml");
    EXPECT_EQ(doc.text("analyses.contacts.cutoff"), "25.0");
    EXPECT_EQ(doc.text("visual.zoom"), "0.01");
    EXPECT_EQ(doc.text("output.video.width"), "7682");
    EXPECT_EQ(doc.text("inputs.dt_ps"), "0.0123456789");
    EXPECT_EQ(doc.text("analyses.energy.target_value"), "-1234567.5");
    EXPECT_EQ(doc.text("inputs.stop"), "");
    EXPECT_EQ(doc.text("visual.base_color"), "tab:orange");
    EXPECT_EQ(doc.text("output.charts.formats"), "jpg, svg");
    EXPECT_EQ(doc.text("visual.renderer"), "vulkan");
    EXPECT_EQ(doc.text("visual.rotation"), "0.0, 0.0, 0.0");
    EXPECT_EQ(doc.text("analyses.use_pbc"), "false");
    EXPECT_EQ(doc.text("layout.panels.3.options"), "{\"mode\": \"instant\", \"transpose\": true}");
    EXPECT_EQ(doc.text("layout.panels.0.options"), "");
    EXPECT_EQ(doc.text("layout.panels.3.options.mode"), "instant");
    EXPECT_EQ(doc.text("layout.panels.3.options.transpose"), "true");
    for (const std::string& path : fieldPaths(doc.project())) {
        if (const auto info = doc.info(path); info->kind != FieldKind::Records && info->kind != FieldKind::Record &&
                                              info->kind != FieldKind::NumberRows) {
            doc.setText(path, doc.text(path));
        }
    }
    EXPECT_FALSE(doc.isModified());
    const FixedClock clock(savedStamp(dir / "p.mdatk.toml"));
    doc.save();
    EXPECT_EQ(readFile(dir / "p.mdatk.toml"), before);
    // the problems flag the unknown values instead of replacing them
    std::vector<std::string> fields;
    for (const Problem& problem : doc.problems()) {
        fields.push_back(problem.field);
    }
    EXPECT_NE(std::find(fields.begin(), fields.end(), "visual.color_mode"), fields.end());
    EXPECT_EQ(std::find(fields.begin(), fields.end(), "visual.renderer"), fields.end());  // RN-20
}

// ---- fields ----------------------------------------------------------------------------------

TEST(ProjectDocument, FieldInfoAndPaths) {
    const Project p;
    EXPECT_EQ(fieldInfo(p, "analyses.contacts.cutoff")->kind, FieldKind::Number);
    EXPECT_EQ(fieldInfo(p, "analyses.contacts.cutoff")->modelPath, "analyses.contacts.cutoff");
    EXPECT_EQ(*fieldInfo(p, "analyses.contacts.cutoff")->comment, "Contact distance cutoff (Å)");
    EXPECT_EQ(fieldInfo(p, "inputs.stop")->kind, FieldKind::OptionalInteger);
    EXPECT_EQ(fieldInfo(p, "inputs.dt_ps")->kind, FieldKind::OptionalNumber);
    EXPECT_EQ(fieldInfo(p, "inputs.trajectories")->kind, FieldKind::TextList);
    EXPECT_EQ(fieldInfo(p, "groups")->kind, FieldKind::Records);
    EXPECT_EQ(fieldInfo(p, "groups.1")->kind, FieldKind::Record);
    EXPECT_EQ(fieldInfo(p, "groups.1.name")->kind, FieldKind::Text);
    EXPECT_EQ(fieldInfo(p, "groups.1.name")->modelPath, "groups.name");
    EXPECT_FALSE(fieldInfo(p, "groups.2.name"));
    EXPECT_FALSE(fieldInfo(p, "groups.x"));
    EXPECT_FALSE(fieldInfo(p, "nope"));
    EXPECT_FALSE(fieldInfo(p, "analyses.contacts.cutoff.x"));
    EXPECT_EQ(fieldInfo(p, "visual.rotation")->kind, FieldKind::NumberList);
    EXPECT_EQ(fieldInfo(p, "density.anchor")->kind, FieldKind::NumberRows);
    EXPECT_EQ(fieldInfo(p, "layout.panels.3.options")->kind, FieldKind::FreeTable);
    EXPECT_TRUE(fieldInfo(p, "layout.panels.3.options")->comment);
    EXPECT_EQ(fieldInfo(p, "layout.panels.3.options.mode")->kind, FieldKind::FreeValue);
    EXPECT_EQ(fieldInfo(p, "layout.panels.3.options.mode")->modelPath, "layout.panels.options.mode");
    EXPECT_FALSE(fieldInfo(p, "layout.panels.3.options.mode")->comment);
    EXPECT_FALSE(fieldInfo(p, "layout.panels.3.options.vmax"));
    EXPECT_EQ(fieldInfo(p, "workflow.parallel")->kind, FieldKind::OptionalInteger);
    EXPECT_EQ(fieldInfo(p, "workflow.supersede_runs")->kind, FieldKind::OptionalBoolean);
    EXPECT_EQ(fieldInfo(p, "")->kind, FieldKind::Record);

    const std::vector<std::string> paths = fieldPaths(p);
    EXPECT_EQ(paths.front(), "notes");
    EXPECT_EQ(paths[1], "inputs.topology");
    EXPECT_NE(std::find(paths.begin(), paths.end(), "groups.1.color"), paths.end());
    EXPECT_NE(std::find(paths.begin(), paths.end(), "layout.panels.3.options"), paths.end());
    EXPECT_EQ(std::find(paths.begin(), paths.end(), "layout.panels.3.options.mode"), paths.end());
    EXPECT_EQ(std::find(paths.begin(), paths.end(), "groups"), paths.end());
    for (const std::string& path : paths) {
        EXPECT_TRUE(fieldInfo(p, path)) << path;
    }
}

TEST(ProjectDocument, GetAndSetByPath) {
    ProjectDocument doc;
    Events events(doc);
    EXPECT_EQ(doc.get("analyses.contacts.cutoff"), Value(4.5));
    EXPECT_EQ(doc.get("inputs.stop"), Value());
    EXPECT_EQ(doc.get("workflow.parallel"), Value());  // a 0.1.x key the project does not have
    EXPECT_EQ(doc.get("groups.1.name"), Value("Target"));
    EXPECT_THROW(doc.get("groups.2.name"), FieldError);

    EXPECT_TRUE(doc.set("analyses.contacts.cutoff", 25));  // an int for a number: 25.0, no clamping
    EXPECT_EQ(doc.project().analyses.contacts.cutoff, 25.0);
    EXPECT_TRUE(doc.isModified());
    EXPECT_EQ(doc.undoLabel(), "Set Analyses › Pairwise contacts › Cutoff");
    ASSERT_EQ(events.seen.size(), 1U);
    EXPECT_EQ(events.seen[0].kind, DocumentChange::Kind::Edit);
    EXPECT_EQ(events.seen[0].fields, std::vector<std::string>{"analyses.contacts.cutoff"});
    EXPECT_TRUE(events.seen[0].modifiedChanged);
    EXPECT_TRUE(events.seen[0].affects("analyses.contacts"));
    EXPECT_TRUE(events.seen[0].affects("analyses.contacts.cutoff"));
    EXPECT_FALSE(events.seen[0].affects("analyses.contacts.group_a"));
    EXPECT_FALSE(events.seen[0].affects("analyses.contact"));
    EXPECT_FALSE(doc.set("analyses.contacts.cutoff", 25.0));  // the same value: nothing
    EXPECT_EQ(events.seen.size(), 1U);

    try {
        doc.set("analyses.contacts.cutoff", "far");
        ADD_FAILURE();
    } catch (const FieldError& error) {
        EXPECT_STREQ(error.what(), "expected a number, got 'far'");
        EXPECT_EQ(error.field(), "analyses.contacts.cutoff");
        EXPECT_EQ(error.problem(&doc.project()).text(), "Analyses › Pairwise contacts › Cutoff: expected a number, got 'far'");
    }
    try {
        doc.set("inputs.step", 2.5);
        ADD_FAILURE();
    } catch (const FieldError& error) {
        EXPECT_STREQ(error.what(), "expected an integer, got 2.5");
    }
    try {
        doc.set("notes", 5);
        ADD_FAILURE();
    } catch (const FieldError& error) {
        EXPECT_STREQ(error.what(), "expected a string, got 5");
    }
    try {
        doc.set("analyses.contact.cutoff", 1.0);
        ADD_FAILURE();
    } catch (const FieldError& error) {
        EXPECT_STREQ(error.what(), "no setting 'analyses.contact.cutoff'");
    }
    EXPECT_THROW(doc.set("groups.2", Table{}), FieldError);
    EXPECT_THROW(doc.set("bogus", 1), FieldError);
    EXPECT_EQ(events.seen.size(), 1U);  // nothing changed by the errors

    // free-form tables: change, add and remove entries
    EXPECT_TRUE(doc.set("layout.panels.3.options.mode", "frame"));
    EXPECT_TRUE(doc.set("layout.panels.3.options.frame", 12));
    EXPECT_EQ(doc.project().layout.panels[3].options,
              (Table{{"mode", "frame"}, {"transpose", true}, {"frame", std::int64_t{12}}}));
    EXPECT_TRUE(doc.set("layout.panels.3.options.transpose", Value()));
    EXPECT_EQ(doc.project().layout.panels[3].options, (Table{{"mode", "frame"}, {"frame", std::int64_t{12}}}));
    EXPECT_THROW(doc.set("layout.panels.3.options", 1), FieldError);
    EXPECT_THROW(doc.set("layout.panels.3.series.0.x", 1), FieldError);
    // 0.1.x keys: set and unset
    EXPECT_TRUE(doc.set("workflow.parallel", 3));
    EXPECT_EQ(doc.project().workflow.parallel, std::optional<std::int64_t>(3));
    EXPECT_TRUE(doc.set("workflow.parallel", Value()));
    EXPECT_FALSE(doc.project().workflow.parallel);
    // planner files must be paths (loadProject refuses the project otherwise)
    doc.edit("Add a step", [](Project& p) { p.workflow.nodes.push_back(PlanNode{}); });
    EXPECT_TRUE(doc.set("workflow.nodes.0.files.structure", "model.pdb"));
    EXPECT_THROW(doc.set("workflow.nodes.0.files.structure", 3), FieldError);
    // the whole project
    EXPECT_TRUE(doc.set("", toPlain(Project{})));
    EXPECT_EQ(doc.project(), Project{});
    EXPECT_THROW(doc.set("", 1), FieldError);
}

// Form errors (AP-10): located messages instead of Python's, the field flagged; nothing changes.
TEST(ProjectDocument, TextOfFields) {
    ProjectDocument doc;
    const auto error = [&](std::string_view path, std::string_view text) {
        try {
            doc.setText(path, text);
        } catch (const FieldError& e) {
            return e.problem(&doc.project()).text();
        }
        return std::string("(no error)");
    };
    EXPECT_EQ(error("inputs.stop", "x"), "Input & cleaning › Frames › Stop: 'x' is not a whole number");
    EXPECT_EQ(error("inputs.stop", "1.5"), "Input & cleaning › Frames › Stop: '1.5' is not a whole number");
    EXPECT_EQ(error("analyses.contacts.cutoff", "4,5"),
              "Analyses › Pairwise contacts › Cutoff: '4,5' is not a number");
    EXPECT_EQ(error("analyses.use_pbc", "maybe"), "Analyses › Use pbc: 'maybe' is not true or false");
    EXPECT_EQ(error("visual.rotation", "1, x, 3"), "Structure view › Rotation: 'x' is not a number");
    EXPECT_EQ(error("layout.panels.3.options", "{bad"),
              "Layout & panels › Panel 4 (Contact map) › Options: not valid JSON: Expecting property name "
              "enclosed in double quotes: line 1 column 2 (char 1)");
    EXPECT_EQ(error("layout.panels.3.options", "[1]"),
              "Layout & panels › Panel 4 (Contact map) › Options: options must be a JSON object");
    EXPECT_FALSE(doc.isModified());

    EXPECT_TRUE(doc.setText("inputs.stop", " 120 "));
    EXPECT_EQ(doc.project().inputs.stop, std::optional<std::int64_t>(120));
    EXPECT_TRUE(doc.setText("inputs.stop", ""));
    EXPECT_FALSE(doc.project().inputs.stop);
    EXPECT_TRUE(doc.setText("inputs.dt_ps", "1e-3"));
    EXPECT_EQ(doc.text("inputs.dt_ps"), "0.001");
    EXPECT_TRUE(doc.setText("analyses.contacts.cutoff", "nan"));
    EXPECT_TRUE(std::isnan(doc.project().analyses.contacts.cutoff));
    EXPECT_FALSE(doc.setText("analyses.contacts.cutoff", "NaN"));  // NaN is NaN: no change
    EXPECT_TRUE(doc.setText("analyses.use_pbc", "True"));
    EXPECT_TRUE(doc.project().analyses.use_pbc);
    EXPECT_TRUE(doc.setText("inputs.trajectories", " a.nc, ,b.nc ,"));
    EXPECT_EQ(doc.project().inputs.trajectories, (std::vector<std::string>{"a.nc", "b.nc"}));
    EXPECT_TRUE(doc.setText("notes", "  kept as typed  "));
    EXPECT_EQ(doc.project().notes, "  kept as typed  ");
    EXPECT_TRUE(doc.setText("layout.panels.3.options", ""));
    EXPECT_EQ(doc.project().layout.panels[3].options, Table{});
    EXPECT_TRUE(doc.setText("layout.panels.3.options", "{\"vmax\": 2, \"mode\": \"frame\"}"));
    EXPECT_EQ(doc.text("layout.panels.3.options.vmax"), "2");
    EXPECT_TRUE(doc.setText("layout.panels.3.options.vmax", "2.5"));
    EXPECT_EQ(doc.get("layout.panels.3.options.vmax"), Value(2.5));
    EXPECT_TRUE(doc.setText("layout.panels.3.options.mode", "instant"));  // a text stays text
    EXPECT_EQ(doc.get("layout.panels.3.options.mode"), Value("instant"));
    EXPECT_THROW(doc.text("groups"), FieldError);
    EXPECT_THROW(doc.text("density.anchor"), FieldError);
}

TEST(ProjectDocument, ListItems) {
    ProjectDocument doc;
    doc.insertItem("groups", 2, Table{{"name", "Ion"}, {"selection", "resname MG"}});
    ASSERT_EQ(doc.project().groups.size(), 3U);
    EXPECT_EQ(doc.project().groups[2].name, "Ion");
    EXPECT_EQ(doc.project().groups[2].color, "#2a78d6");  // defaults for what the item leaves out
    doc.moveItem("groups", 2, 0);
    EXPECT_EQ(groupNames(doc.project()), (std::vector<std::string>{"Ion", "Aptamer", "Target"}));
    doc.removeItem("groups", 1);
    EXPECT_EQ(groupNames(doc.project()), (std::vector<std::string>{"Ion", "Target"}));
    doc.insertItem("inputs.trajectories", 0, "a.nc");
    EXPECT_THROW(doc.insertItem("inputs.trajectories", 0, 1), FieldError);
    EXPECT_THROW(doc.insertItem("inputs.trajectories", 5, "b.nc"), FieldError);
    EXPECT_THROW(doc.removeItem("notes", 0), FieldError);
    EXPECT_THROW(doc.moveItem("groups", 0, 2), FieldError);
    EXPECT_EQ(doc.undoLabel(), "Add to Input & cleaning › Input files › Trajectories");
    doc.undo();
    doc.undo();
    doc.undo();
    doc.undo();
    EXPECT_EQ(doc.project(), Project{});
    EXPECT_FALSE(doc.isModified());
}

// ---- undo ------------------------------------------------------------------------------------

TEST(ProjectDocument, UndoAndRedo) {
    TempDir dir;
    ProjectDocument doc;
    Events events(doc);
    EXPECT_FALSE(doc.canUndo());
    doc.setText("analyses.contacts.cutoff", "5");
    doc.setText("analyses.contacts.cutoff", "6");
    doc.setText("inputs.topology", "/data/a.parm7", "Choose topology");
    EXPECT_EQ(doc.undoLabel(), "Choose topology");
    doc.saveAs(dir / "p.mdatk.toml");
    EXPECT_FALSE(doc.isModified());
    doc.undo();
    EXPECT_TRUE(doc.isModified());
    EXPECT_EQ(doc.project().inputs.topology, "");
    EXPECT_EQ(doc.redoLabel(), "Choose topology");
    EXPECT_EQ(events.seen.back().kind, DocumentChange::Kind::Undo);
    EXPECT_EQ(events.seen.back().fields, std::vector<std::string>{"inputs.topology"});
    EXPECT_TRUE(events.seen.back().modifiedChanged);
    doc.undo();
    EXPECT_EQ(doc.project().analyses.contacts.cutoff, 5.0);
    doc.redo();
    doc.redo();
    EXPECT_FALSE(doc.isModified());  // back at the saved state
    EXPECT_EQ(doc.project().inputs.topology, "/data/a.parm7");
    EXPECT_FALSE(doc.canRedo());
    doc.undo();
    doc.setText("notes", "x");  // a new change drops the redo steps
    EXPECT_FALSE(doc.canRedo());
    EXPECT_TRUE(doc.isModified());
    doc.undo();
    doc.undo();
    doc.undo();
    EXPECT_FALSE(doc.canUndo());
    EXPECT_EQ(doc.project().analyses.contacts.cutoff, 4.5);
    doc.undo();  // nothing left: nothing happens
    EXPECT_EQ(doc.project().analyses.contacts.cutoff, 4.5);
}

// A drag is one step (D5): the changes inside a group are applied and notified at once, undone
// together; cancelGroup restores the state before the drag.
TEST(ProjectDocument, GroupsAreOneStep) {
    ProjectDocument doc;
    Events events(doc);
    doc.beginGroup("Rotate structure");
    for (int i = 1; i <= 10; ++i) {
        doc.set("visual.rotation", Array{Value(i * 1.0), Value(0.0), Value(0.0)});
        doc.beginGroup("inner");  // nested: the outer label counts
        doc.set("visual.zoom", 1.0 + i * 0.1);
        doc.endGroup();
    }
    EXPECT_TRUE(doc.inGroup());
    doc.endGroup();
    EXPECT_FALSE(doc.inGroup());
    EXPECT_EQ(events.seen.size(), 20U);
    EXPECT_EQ(events.seen[0].fields, std::vector<std::string>{"visual.rotation.0"});
    EXPECT_EQ(doc.undoLabel(), "Rotate structure");
    doc.undo();
    EXPECT_FALSE(doc.canUndo());
    EXPECT_EQ(doc.project(), Project{});
    EXPECT_FALSE(doc.isModified());
    doc.redo();
    EXPECT_EQ(doc.project().visual.rotation, (std::vector<double>{10.0, 0.0, 0.0}));

    doc.beginGroup("Move map");
    doc.set("density.fit_info", "moved");
    doc.cancelGroup();
    EXPECT_EQ(doc.project().density.fit_info, "");
    EXPECT_EQ(doc.undoLabel(), "Rotate structure");
    EXPECT_FALSE(doc.inGroup());

    doc.beginGroup("Nothing");
    doc.endGroup();
    EXPECT_EQ(doc.undoLabel(), "Rotate structure");
    doc.beginGroup("Typed");  // undo during a group ends it first
    doc.setText("notes", "abc");
    doc.undo();
    EXPECT_EQ(doc.project().notes, "");
    EXPECT_EQ(doc.redoLabel(), "Typed");
}

TEST(ProjectDocument, UndoKeepsAtMostTheLimit) {
    ProjectDocument doc;
    for (std::size_t i = 0; i < ProjectDocument::kUndoLimit + 5; ++i) {
        doc.set("inputs.start", static_cast<std::int64_t>(i + 1));
    }
    std::size_t steps = 0;
    while (doc.canUndo()) {
        doc.undo();
        ++steps;
    }
    EXPECT_EQ(steps, ProjectDocument::kUndoLimit);
    EXPECT_EQ(doc.project().inputs.start, 5);
    EXPECT_TRUE(doc.isModified());  // the saved state was dropped
}

TEST(ProjectDocument, EditAndListeners) {
    ProjectDocument doc;
    std::vector<std::string> order;
    const auto first = doc.addListener([&](const DocumentChange&) { order.push_back("first"); });
    ProjectDocument::ListenerId second = 0;
    second = doc.addListener([&](const DocumentChange&) {
        order.push_back("second");
        doc.removeListener(second);  // from the next event on
    });
    EXPECT_FALSE(doc.edit("Nothing", [](Project&) {}));
    EXPECT_TRUE(doc.edit("Swap panels", [](Project& p) { std::swap(p.layout.panels[1], p.layout.panels[2]); }));
    EXPECT_EQ(order, (std::vector<std::string>{"first", "second"}));
    doc.undo();
    EXPECT_EQ(order, (std::vector<std::string>{"first", "second", "first"}));
    doc.removeListener(first);
    EXPECT_THROW(doc.edit("Fails", [](Project& p) {
        p.notes = "changed";
        throw std::runtime_error("no");
    }),
                 std::runtime_error);
    EXPECT_EQ(doc.project().notes, "");
    EXPECT_EQ(order.size(), 3U);

    Events events(doc);
    doc.newProject();
    ASSERT_EQ(events.seen.size(), 1U);
    EXPECT_EQ(events.seen[0].kind, DocumentChange::Kind::Replaced);
    EXPECT_TRUE(events.seen[0].affects("anything.at.all"));
    EXPECT_FALSE(doc.canUndo());
    EXPECT_FALSE(doc.canRedo());
}

// ---- files -----------------------------------------------------------------------------------

// D8: Save As to another folder rebases relative paths so they name the same files; absolute and
// "~" paths stay; the undo steps follow.
TEST(ProjectDocument, SaveAsRebasesRelativePaths) {
    TempDir dir;
    fs::create_directories(dir / "a");
    fs::create_directories(dir / "b" / "c");
    Project p;
    p.inputs.topology = "inputs/sys.parm7";
    p.inputs.trajectories = {"../shared/run.nc", "/abs/run2.nc"};
    p.density.map_path = "~/maps/m.mrc";
    p.output.out_dir = "MDATK_output";
    PlanNode node;
    node.tool = "input-file";
    node.files = Table{{"file", "a.pdb; ../b.pdb;/x/c.pdb"}};
    node.settings = Table{{"path", "inputs/model.pdb"}, {"kind", "pdb"}};
    p.workflow.nodes.push_back(node);
    saveProject(dir / "a" / "p.mdatk.toml", p);

    ProjectDocument doc;
    doc.open(dir / "a" / "p.mdatk.toml");
    const Project before = doc.resolvedProject();
    doc.setText("notes", "edited");
    Events events(doc);
    doc.saveAs(dir / "a" / "q.mdatk.toml");  // the same folder: nothing to rebase
    EXPECT_EQ(doc.project().inputs.topology, "inputs/sys.parm7");
    EXPECT_TRUE(events.seen.back().fields.empty());
    doc.saveAs(dir / "b" / "c" / "p.mdatk.toml");
    const Project& q = doc.project();
    EXPECT_EQ(q.inputs.topology, "../../a/inputs/sys.parm7");
    EXPECT_EQ(q.inputs.trajectories, (std::vector<std::string>{"../../shared/run.nc", "/abs/run2.nc"}));
    EXPECT_EQ(q.density.map_path, "~/maps/m.mrc");
    EXPECT_EQ(q.output.out_dir, "../../a/MDATK_output");
    EXPECT_EQ(q.workflow.nodes[0].files, (Table{{"file", "../../a/a.pdb;../../b.pdb;/x/c.pdb"}}));
    EXPECT_EQ(q.workflow.nodes[0].settings, (Table{{"path", "../../a/inputs/model.pdb"}, {"kind", "pdb"}}));
    EXPECT_EQ(events.seen.back().kind, DocumentChange::Kind::Saved);
    EXPECT_EQ(events.seen.back().fields,
              (std::vector<std::string>{"inputs.topology", "inputs.trajectories.0", "output.out_dir",
                                        "workflow.nodes.0.files.file", "workflow.nodes.0.settings.path"}));
    // the same files as before (the planner files as resolvePaths joins them)
    const Project after = doc.resolvedProject();
    EXPECT_EQ(after.inputs, before.inputs);
    EXPECT_EQ(after.output.out_dir, before.output.out_dir);
    EXPECT_EQ(after.workflow.nodes[0].files, before.workflow.nodes[0].files);
    EXPECT_EQ(loadProject(dir / "b" / "c" / "p.mdatk.toml").inputs, before.inputs);
    doc.undo();  // the undo step was rebased too
    EXPECT_EQ(doc.project().inputs.topology, "../../a/inputs/sys.parm7");
    EXPECT_EQ(doc.project().notes, "");
}

TEST(ProjectDocument, RelativePath) {
#ifndef _WIN32
    EXPECT_EQ(relativePath("/a/b/c", "/a/b"), "c");
    EXPECT_EQ(relativePath("/a/b", "/a/b"), ".");
    EXPECT_EQ(relativePath("/a/x/y", "/a/b/c"), "../../x/y");
    EXPECT_EQ(relativePath("/", "/a/b"), "../..");
    EXPECT_EQ(relativePath("/a//b/./c", "/a/b/"), "c");
#else
    EXPECT_EQ(relativePath("C:/a/b/c", "C:/a/b"), "c");
    EXPECT_EQ(relativePath("c:/a/x", "C:/a/b"), "../x");
    EXPECT_EQ(relativePath("C:/Work/Proj/inputs/a.nc", "c:/work/proj"), "inputs/a.nc");  // names ignore case
    EXPECT_EQ(relativePath("C:/work/Other/b", "C:/Work/proj"), "../Other/b");
    EXPECT_EQ(relativePath("D:/a", "C:/a"), std::nullopt);
    EXPECT_EQ(relativePath("//srv/share/x", "//srv/share/y"), "../x");
    EXPECT_EQ(relativePath("//srv/other/x", "//srv/share/y"), std::nullopt);
#endif
}

// The step log a run appends is no change on disk; an edit of the settings by another program is.
TEST(ProjectDocument, ChangedOnDisk) {
    TempDir dir;
    const fs::path path = dir / "p.mdatk.toml";
    saveProject(path, Project{});
    ProjectDocument doc;
    EXPECT_FALSE(doc.changedOnDisk());
    doc.open(path);
    EXPECT_FALSE(doc.changedOnDisk());
    appendLog(path, LogEntry("render", "Render"));
    EXPECT_FALSE(doc.changedOnDisk());
    Project other;
    other.notes = "from MDATK";
    saveProject(path, other);
    EXPECT_TRUE(doc.changedOnDisk());
    doc.save();
    EXPECT_FALSE(doc.changedOnDisk());
    fs::remove(path);
    EXPECT_TRUE(doc.changedOnDisk());
}

TEST(ProjectDocument, OpenErrorsLeaveTheDocument) {
    ProjectDocument doc;
    doc.setText("notes", "mine");
    EXPECT_THROW(doc.open(kProjects / "broken.mdatk.toml"), ProjectError);
    EXPECT_THROW(doc.open(kProjects / "other_format.toml"), ProjectError);
    EXPECT_THROW(doc.open(kProjects / "nope.mdatk.toml"), ProjectError);
    EXPECT_EQ(doc.project().notes, "mine");
    EXPECT_TRUE(doc.canUndo());
    EXPECT_THROW(doc.save(), ProjectError);  // no file yet
    std::vector<std::string> warnings;
    TempDir dir;
    writeFile(dir / "u.mdatk.toml", "format = \"mdatk-project\"\nbogus = 1\n[inputs]\nnope = 2\n");
    doc.open(dir / "u.mdatk.toml", &warnings);
    EXPECT_EQ(warnings.size(), 1U);  // D9: one warning on open
    EXPECT_FALSE(doc.canUndo());
    EXPECT_TRUE(doc.setText("notes", "x"));
    doc.save();
    EXPECT_NE(readFile(dir / "u.mdatk.toml").find("bogus = 1"), std::string::npos);  // kept (PRJ-09)
}

// ---- residue groups (MDATK 0.1.1) -----------------------------------------------------------

TEST(ProjectDocument, RenameGroupAsMdatkLatest) {
    const Value cases = readJson(kCases);
    std::size_t count = 0;
    for (const Value& c : items(cases, "rename")) {
        const std::string name = c.find("name")->stringOr();
        Project p = projectFromPlain(tagged(*c.find("input")));
        const std::string oldName = c.find("old")->stringOr();
        const std::string newName = c.find("new")->stringOr();
        if (const Value* error = c.find("error")) {
            const std::string type = error->find("type")->stringOr();
            try {
                renameGroup(p, oldName, newName);
                ADD_FAILURE() << name;
            } catch (const std::out_of_range& e) {
                EXPECT_EQ(type, "KeyError") << name;
                EXPECT_EQ(e.what(), error->find("message")->stringOr()) << name;
            } catch (const ValueError& e) {
                EXPECT_EQ(type, "ValueError") << name;
                EXPECT_EQ(e.what(), error->find("message")->stringOr()) << name;
            }
        } else {
            EXPECT_EQ(renameGroup(p, oldName, newName), dotted(strings(*c.find("changed")))) << name;
            EXPECT_EQ(pythonRepr(toPlain(p)), pythonRepr(tagged(*c.find("output")))) << name;
        }
        ++count;
    }
    EXPECT_EQ(count, 14U);
    for (const Value& c : items(cases, "references")) {
        const Project p = projectFromPlain(tagged(*c.find("input")));
        EXPECT_EQ(groupReferences(p, c.find("group")->stringOr()), dotted(strings(*c.find("found"))))
            << c.find("name")->stringOr();
    }
}

TEST(ProjectDocument, RenameGroupIsOneStep) {
    ProjectDocument doc;
    Events events(doc);
    const std::vector<std::string> changed = doc.renameGroup("Aptamer", "Apt");
    EXPECT_EQ(changed.size(), 9U);
    EXPECT_EQ(doc.project().groups[0].name, "Apt");
    EXPECT_EQ(doc.project().analyses.contacts.group_a, "Apt");
    EXPECT_EQ(doc.undoLabel(), "Rename group 'Aptamer' to 'Apt'");
    ASSERT_EQ(events.seen.size(), 1U);
    EXPECT_TRUE(events.seen[0].affects("groups.0.name"));
    EXPECT_TRUE(events.seen[0].affects("analyses.rmsd.targets.0.selection"));
    EXPECT_EQ(doc.groupReferences("Apt").size(), 9U);
    EXPECT_TRUE(doc.renameGroup("Apt", "Apt").empty());
    EXPECT_EQ(events.seen.size(), 1U);
    EXPECT_THROW(doc.renameGroup("Apt", "Target"), ValueError);
    EXPECT_THROW(doc.renameGroup("Nope", "X"), std::out_of_range);
    doc.undo();
    EXPECT_EQ(doc.project(), Project{});
    EXPECT_EQ(groupNameProblem("ok_1.2-3"), "");
    EXPECT_EQ(groupNameProblem(" "), "A residue group needs a name.");
    EXPECT_EQ(groupNameProblem("Target", {"Target"}), "There is already a residue group named 'Target'.");
}

// ---- problems --------------------------------------------------------------------------------

// Every validate() problem is located on its field; where MDATK 0.1.1 has the same check, field,
// location and message are 0.1.1's (D7: which problems there are stays validate()'s).
TEST(ProjectDocument, ProblemsOfEveryValidateCaseAreLocated) {
    const Value vectors = readJson(kData / "project-model" / "validate_cases.json");
    const Value latest = readJson(kCases);
    TempDir dir;
    for (const std::string& name : strings(*vectors.find("files"))) {
        writeFile(dir / name, "");
    }
    for (const std::string& name : strings(*vectors.find("folders"))) {
        fs::create_directory(dir / name);
    }
    const ScopedCwd cwd(dir.path());
    const Array& expected = items(latest, "problems");
    std::size_t index = 0;
    std::size_t matched = 0;
    std::size_t located = 0;
    for (const Value& c : items(vectors, "cases")) {
        const std::string name = c.find("name")->stringOr();
        const Project p = projectFromPlain(tagged(*c.find("input")));
        const std::vector<Problem> problems = projectProblems(p);
        std::vector<std::string> checks;
        for (const Problem& problem : problems) {
            if (problem.level == "error") {
                checks.push_back(problem.check);
                EXPECT_FALSE(problem.field.empty()) << name << ": " << problem.check;
                EXPECT_TRUE(fieldInfo(p, problem.field)) << name << ": " << problem.field;
                EXPECT_FALSE(problem.page().empty()) << name << ": " << problem.field;
                EXPECT_EQ(problem.location, describeField(problem.field, &p)) << name;
                EXPECT_EQ(problem.text(), problem.location + ": " + problem.message);
                ++located;
            } else {
                EXPECT_EQ(problem.field, "layout.panels");
            }
        }
        EXPECT_EQ(checks, validate(p)) << name;
        const Value& mdatk = expected[index++];
        ASSERT_EQ(mdatk.find("name")->stringOr(), name);
        if (const Value* list = mdatk.find("problems")) {
            for (const Value& m : *list->get<Array>()) {
                const std::string message = m.find("message")->stringOr();
                const std::string field = m.find("field")->stringOr();
                const bool known = std::any_of(problems.begin(), problems.end(), [&](const Problem& x) {
                    return x.message == message && x.level == m.find("level")->stringOr();
                });
                if (!known) {
                    continue;  // a check MDBridge does not have, or words it changed
                }
                const auto same = std::find_if(problems.begin(), problems.end(), [&](const Problem& x) {
                    return x.message == message && x.field == field;
                });
                ASSERT_NE(same, problems.end()) << name << ": " << field << ": " << message;
                EXPECT_EQ(same->location, m.find("location")->stringOr()) << name << ": " << message;
                ++matched;
            }
        }
    }
    EXPECT_EQ(index, expected.size());
    EXPECT_GT(located, 150U);
    EXPECT_GT(matched, 100U);
}

TEST(ProjectDocument, DescribeAsMdatkLatest) {
    const Value cases = readJson(kCases);
    const Project p = projectFromPlain(tagged(*items(cases, "rename")[0].find("input")));
    for (const Value& c : items(cases, "describe")) {
        const std::string field = c.find("field")->stringOr();
        EXPECT_EQ(describeField(field), c.find("plain")->stringOr()) << field;
        EXPECT_EQ(describeField(field, &p), c.find("project")->stringOr()) << field;
    }
    EXPECT_EQ(problemPage("analyses.contacts.cutoff"), "analyses");
    EXPECT_EQ(problemPage("density.map_path"), "structure");
    EXPECT_EQ(problemPage("workflow.nodes"), "");
    EXPECT_EQ(problemPageTitle("outputs"), "Outputs");
}

TEST(ProjectDocument, DocumentProblemsUseResolvedPaths) {
    TempDir dir;
    writeFile(dir / "sys.parm7", "");
    writeFile(dir / "run.nc", "");
    Project p;
    p.inputs.topology = "sys.parm7";
    p.inputs.trajectories = {"run.nc", "gone.nc"};
    saveProject(dir / "p.mdatk.toml", p);
    ProjectDocument doc;
    doc.open(dir / "p.mdatk.toml");
    const std::vector<Problem> problems = doc.problems();
    ASSERT_EQ(problems.size(), 1U);
    EXPECT_EQ(problems[0].field, "inputs.trajectories");
    EXPECT_EQ(problems[0].text(), "Input & cleaning › Input files › Trajectories: Trajectory file not found: " +
                                      resolvePath("gone.nc", doc.folder()));
}

// ---- recovery copies (D6) --------------------------------------------------------------------

TEST(ProjectDocument, RecoveryCopies) {
    TempDir dir;
    const fs::path recovery = dir / "recovery";
    fs::create_directories(dir / "proj");
    Project p;
    p.inputs.topology = "sys.parm7";
    saveProject(dir / "proj" / "p.mdatk.toml", p);
    const std::string saved = readFile(dir / "proj" / "p.mdatk.toml");

    ProjectDocument doc;
    doc.open(dir / "proj" / "p.mdatk.toml");
    doc.setText("notes", "unsaved work");
    const fs::path file = writeRecovery(doc, recovery);
    EXPECT_EQ(file, recoveryPath(doc, recovery));
    EXPECT_EQ(readFile(dir / "proj" / "p.mdatk.toml"), saved);  // the project file is never written
    ProjectDocument untitled;
    untitled.setText("inputs.topology", "rel.parm7");
    const fs::path second = writeRecovery(untitled, recovery);
    const std::string untitledTopology = untitled.resolvedProject().inputs.topology;
    EXPECT_NE(second, file);
    writeFile(recovery / "recovery-junk.mdatk.toml", "not toml [");
    writeFile(recovery / "other.txt", "");

    const std::vector<RecoveryFile> found = findRecoveries(recovery);
    ASSERT_EQ(found.size(), 2U);
    const auto mine = std::find_if(found.begin(), found.end(), [&](const RecoveryFile& r) { return r.file == file; });
    ASSERT_NE(mine, found.end());
    EXPECT_EQ(mine->original, doc.path());
    EXPECT_EQ(mine->folder, doc.folder());

    ProjectDocument restored;
    openRecovery(restored, *mine);
    EXPECT_EQ(restored.project(), doc.project());
    EXPECT_EQ(restored.path(), doc.path());
    EXPECT_TRUE(restored.isModified());
    EXPECT_FALSE(restored.canUndo());
    restored.save();
    EXPECT_FALSE(restored.isModified());
    EXPECT_TRUE(removeRecovery(file));
    EXPECT_FALSE(removeRecovery(file));
    EXPECT_EQ(findRecoveries(recovery).size(), 1U);

    // a recovery of a new project made in another folder: its relative paths are rebased
    const auto other = findRecoveries(recovery).front();
    EXPECT_TRUE(other.original.empty());
    fs::create_directories(dir / "elsewhere");
    const ScopedCwd cwd(dir / "elsewhere");
    ProjectDocument fresh;
    openRecovery(fresh, other);
    EXPECT_FALSE(fresh.hasPath());
    EXPECT_EQ(fresh.resolvedProject().inputs.topology, untitledTopology);
    EXPECT_NE(fresh.project().inputs.topology, "rel.parm7");
    EXPECT_TRUE(findRecoveries(dir / "missing").empty());
}
