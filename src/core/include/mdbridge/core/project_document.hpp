#pragma once
// The edit model of the GUI (milestone 5, W0-B): one open project, edited by path, with undo, a
// modified flag, change notifications and located problems. No Qt; the GUI's field widgets bind to
// ProjectDocument by field path, and so do tests.
//
// Field paths are the project-file keys joined with '.', list items by their 0-based index:
// "analyses.contacts.cutoff", "groups.1.name", "inputs.trajectories.0", "layout.panels.2.options"
// and, inside the free-form tables (panel options, planner settings and files), their keys at any
// depth: "layout.panels.2.options.mode". "" is the whole project. These are the paths of MDATK's
// GUI bindings (gui/widgets.py get_path / set_path) and of MDATK 0.1.1's problems (problems.py).
//
// Lossless (milestone 5 AP-01 .. AP-05): the document holds exactly what was read. A value is
// never clamped to a range, rounded, re-spelt or replaced by a default; text() of a number is its
// shortest round-trip text (the project writer's), colours and unknown choices are text kept as
// they are, unknown keys stay where they were (D9, PRJ-09). Setting a field to the text it already
// shows changes nothing, so a field the user does not edit is saved back as it was read.
//
// Paths as written (AP-06, owner decision D8): the document keeps file paths as the file has them
// (relative to the project's folder, "~", absolute); resolvedProject() is the project with them
// made absolute, as loadProject() gives it, for running. Save keeps them; Save As to another folder
// rebases the relative ones so that they still name the same files.
//
// Undo (owner decision D5): every change is one step holding a snapshot of the whole project
// (Project values are small); beginGroup / endGroup make a drag one step. Undo is not saved.
//
// Not thread-safe: use one document from one thread (the GUI's); workers get copies
// (resolvedProject()).

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// ---- fields ----------------------------------------------------------------------------------

// What a field path holds, from the model (not from the current value).
enum class FieldKind {
    Text,             // std::string
    Integer,          // std::int64_t
    Number,           // double
    OptionalInteger,  // std::optional<std::int64_t> (unset: Null, text "")
    OptionalNumber,   // std::optional<double>
    OptionalText,     // std::optional<std::string> (PlanNode.adopted, folder_tag)
    OptionalBoolean,  // std::optional<bool> (WorkflowParams.supersede_runs)
    Boolean,          // bool
    TextList,         // std::vector<std::string> (trajectories, groups, series, formats ...)
    NumberList,       // std::vector<double> (visual.rotation, density.transform)
    NumberRows,       // std::vector<std::vector<double>> (density.anchor)
    Records,          // a list of tables of the model (groups, layout.panels, ...)
    Record,           // a table of the model (analyses.contacts, one group, one panel)
    FreeTable,        // a free-form table (panel options, planner settings / files) or a table in one
    FreeValue,        // any other value inside a free-form table
};

struct FieldInfo {
    FieldKind kind = FieldKind::Text;
    // The path with list indices left out ("groups.name", "layout.panels.options.mode"), the key
    // of modelComment().
    std::string modelPath;
    // The field's comment in project files (nullopt inside free-form tables).
    std::optional<std::string_view> comment;
};

// The kind of `path` in `project`: nullopt when the path names nothing (no such field, a list
// index past the end, a key a free-form table does not have, an unknown key). A key a free-form
// table does not have yet can still be set().
std::optional<FieldInfo> fieldInfo(const Project& project, std::string_view path);

// Every field of the project down to the values a widget edits, in file order: scalars,
// optionals, lists of scalars and number rows, free-form tables (as one field each); the items of
// lists of records by index ("groups.0.name", "groups.0.residues", ..., "layout.panels.0.options").
std::vector<std::string> fieldPaths(const Project& project);

// A field path that cannot be read or set, or a value or text that does not fit it. what() is the
// message without the location; field() the path; problem(project) the located problem
// ("Input & cleaning › Frames › Stop: 'x' is not a whole number").
class FieldError : public ValueError {
public:
    FieldError(std::string field, const std::string& message) : ValueError(message), m_field(std::move(field)) {}
    const std::string& field() const { return m_field; }
    Problem problem(const Project* project = nullptr) const;

private:
    std::string m_field;
};

// The text a field shows (Text, Integer, Number, Optional*, Boolean, TextList, NumberList,
// FreeTable, FreeValue): numbers as the project writer writes them (shortest round trip: 4.5,
// 25.0, 1e-05, nan, inf), unset optionals as "", booleans "true" / "false", lists of text as
// "a, b" (MDATK's comma-separated fields), lists of numbers as "0.0, 90.0, 0.0", free-form tables
// as JSON ("" for an empty one: MDATK's slot options field), free values as JSON. Throws FieldError
// for Records, Record and NumberRows ("<path> has no text form") and for paths that name nothing.
std::string fieldText(const Project& project, std::string_view path);

// The value `text` means for `path` (what setText would store), or FieldError with the message
// of a form error (AP-10): "'x' is not a whole number", "'x' is not a number", "'maybe' is not
// true or false", "not valid JSON: Expecting value: line 1 column 1 (char 0)", "options must be a
// JSON object". Surrounding whitespace is ignored for numbers; "" is unset for optionals and the
// empty list / table for lists and tables; text fields take the text as it is.
Value parseFieldText(const Project& project, std::string_view path, std::string_view text);

// ---- changes ---------------------------------------------------------------------------------

struct DocumentChange {
    enum class Kind {
        Edit,     // set / setText / edit / renameGroup (also each change inside a group)
        Undo,
        Redo,
        Replaced, // a new, opened or recovered project: every field may have changed
        Saved,    // save / saveAs: values unchanged unless Save As rebased paths (then listed)
    };
    Kind kind = Kind::Edit;
    // The fields whose values changed (paths as above; a list that grew or shrank is named as the
    // list). {""} for Replaced.
    std::vector<std::string> fields;
    // The step's label ("Set Contact cutoff", "Rotate structure"; "" for Replaced / Saved).
    std::string label;
    // isModified() changed with this event.
    bool modifiedChanged = false;

    // True when the value at `path` may have changed: a changed field is `path`, inside it, or
    // holds it.
    bool affects(std::string_view path) const;
};

// ---- the document ----------------------------------------------------------------------------

class ProjectDocument {
public:
    // The number of this document within the process (1, 2, ...).
    std::uint64_t serial() const { return m_serial; }

    using Listener = std::function<void(const DocumentChange&)>;
    using ListenerId = std::uint64_t;

    // At most this many undo steps are kept (the oldest are dropped).
    static constexpr std::size_t kUndoLimit = 200;

    // A new project: the defaults, no file, not modified.
    ProjectDocument();
    ProjectDocument(const ProjectDocument&) = delete;
    ProjectDocument& operator=(const ProjectDocument&) = delete;

    // ---- the project and its file ----

    // The project as it will be saved (paths as written).
    const Project& project() const { return m_project; }
    // The project with its paths made absolute against folder(), as loadProject() reads the file:
    // what analyses and renders use. Throws ProjectError where loadProject would (a planner file
    // that is no path).
    Project resolvedProject() const;
    // The project file ("" for a new project).
    const std::filesystem::path& path() const { return m_path; }
    bool hasPath() const { return !m_path.empty(); }
    // The folder relative paths are relative to: the project file's (absolute), or the current
    // directory for a new project.
    std::string folder() const;
    // The file is a run's parameter log (opened to replay it): save() refuses, use saveAs().
    bool isParamsLog() const { return m_kind == FileKind::ParamsLog; }

    // Replaces the document (Replaced event; undo cleared; not modified).
    void newProject();
    // Opens a project file or parameter log. Paths stay as written. `warnings` receives the
    // unknown-keys warning (D9: one, on open). Throws ProjectError as loadProject does; the
    // document is unchanged then.
    void open(const std::filesystem::path& path, std::vector<std::string>* warnings = nullptr);
    // Replaces the document with `project` (paths as written) belonging to `path` ("" for none).
    // `modified`: the state counts as unsaved (a recovered copy).
    void replace(Project project, std::filesystem::path path = {}, bool modified = false);

    // Saves to path() (saveProject: the step log is kept). Throws ProjectError without a path or
    // for a parameter log, std::filesystem::filesystem_error when the file cannot be written.
    void save();
    // Saves to `path` and makes it the document's file. Relative paths are rebased from folder()
    // to the new file's folder (rebasePaths), in the undo steps too; with the same folder nothing
    // changes. The new file gets the old file's step log only if it is the same file.
    void saveAs(const std::filesystem::path& path);

    // True when the file changed on disk since it was opened or saved here (size, time or bytes
    // differ; also when it was deleted). False for a new project.
    bool changedOnDisk() const;

    // ---- reading and editing fields ----

    // The plain value at `path` (Null for an unset optional, also for a 0.1.x key the file does not
    // have). Throws FieldError for a path that names nothing.
    Value get(std::string_view path) const;
    std::string text(std::string_view path) const { return fieldText(m_project, path); }
    std::optional<FieldInfo> info(std::string_view path) const { return fieldInfo(m_project, path); }

    // Sets the field at `path` (a model field, a list item, or an entry of a free-form table, which
    // is added when missing; Null removes a free-form entry). The value is checked as a project file
    // entry is (MDATK's messages: "expected a number, got 'far'"; text fields need a string).
    // Returns false and changes nothing when the value equals the current one (NaN equals NaN);
    // else one undo step labelled `label` (default "Set <describeField(path)>"). Throws FieldError;
    // the document is unchanged then.
    bool set(std::string_view path, const Value& value, std::string_view label = {});
    // set(path, parseFieldText(...)), but text equal to text(path) changes nothing (so a value is
    // kept exactly as read until the user types something else). Throws FieldError (AP-10).
    bool setText(std::string_view path, std::string_view text, std::string_view label = {});

    // A structural change (add / remove / move rows, swap panels, preset change, canvas actions):
    // `change` edits the project in place; one undo step unless nothing changed. Listeners get
    // the changed fields. An exception from `change` leaves the document unchanged.
    bool edit(std::string_view label, const std::function<void(Project&)>& change);

    // Lists (Records, TextList, NumberList, NumberRows): insert `item` at `index` (== size appends),
    // remove item `index`, move item `from` to position `to`. Checked as set(). Throws FieldError
    // for a path that is no list or an index out of range.
    void insertItem(std::string_view listPath, std::size_t index, const Value& item, std::string_view label = {});
    void removeItem(std::string_view listPath, std::size_t index, std::string_view label = {});
    void moveItem(std::string_view listPath, std::size_t from, std::size_t to, std::string_view label = {});

    // ---- residue groups (MDATK 0.0.10) ----

    // Project.rename_group of MDATK 0.1.1 as one undo step (see core::renameGroup); returns the
    // fields that changed besides the group's name. Throws std::out_of_range (no such group) or
    // ValueError (empty, invalid or taken name; several groups named `oldName`).
    std::vector<std::string> renameGroup(std::string_view oldName, std::string_view newName);
    // The fields that refer to group `name` (what still uses a group about to be deleted).
    std::vector<std::string> groupReferences(std::string_view name) const;

    // ---- undo ----

    bool canUndo() const { return !m_undo.empty(); }
    bool canRedo() const { return !m_redo.empty(); }
    // Labels of the next undo / redo step ("" if none).
    std::string undoLabel() const;
    std::string redoLabel() const;
    void undo();
    void redo();
    void clearUndo();

    // Groups the changes until the matching endGroup into one undo step labelled `label` (a drag,
    // a typed burst). Changes inside are applied and notified at once. Groups nest (the outermost
    // label counts). An empty group adds no step.
    void beginGroup(std::string_view label);
    void endGroup();
    // Ends the group and restores the project as it was at beginGroup (Esc during a drag).
    void cancelGroup();
    bool inGroup() const { return m_groupDepth > 0; }

    // ---- state ----

    // The project differs from the file's last opened or saved state (undoing back to it clears
    // the flag; a new project is unmodified until edited).
    bool isModified() const { return m_currentId != m_savedId; }

    // Every problem of the project (projectProblems of resolvedProject(); for a project whose paths
    // cannot be resolved, of the project as written).
    std::vector<Problem> problems() const;

    // Listeners are called after each change, in the order they were added. A listener may add or
    // remove listeners (taking effect from the next event) but must not edit the document.
    ListenerId addListener(Listener listener);
    void removeListener(ListenerId id);

private:
    struct Step {
        Project project;
        std::uint64_t id = 0;
        std::string label;
    };

    void commit(Project next, std::string_view label);
    void notify(DocumentChange change);
    void rememberDisk();

    Project m_project;
    std::uint64_t m_currentId = 0;
    std::uint64_t m_savedId = 0;
    std::uint64_t m_nextId = 1;
    std::filesystem::path m_path;
    FileKind m_kind = FileKind::Project;
    std::vector<Step> m_undo;
    std::vector<Step> m_redo;
    int m_groupDepth = 0;
    std::optional<Step> m_groupStart;  // the state at the outermost beginGroup
    std::string m_groupLabel;
    bool m_groupChanged = false;
    std::vector<std::pair<ListenerId, Listener>> m_listeners;
    ListenerId m_nextListener = 1;
    std::uint64_t m_serial = 0;  // number of this document within the process (recovery names)
    // what the file was when last read or written here
    std::optional<std::uintmax_t> m_diskSize;
    std::filesystem::file_time_type m_diskTime{};
    std::string m_diskBytes;
};

// ---- residue groups --------------------------------------------------------------------------

// group_name_problem of MDATK 0.1.1: why `name` cannot name a residue group ("" when it can):
// "A residue group needs a name.", "Group name 'a b' may only contain letters, digits, '_', '.' and
// '-' (no spaces).", "There is already a residue group named 'X'." (`taken`: the other names).
std::string groupNameProblem(std::string_view name, const std::vector<std::string>& taken = {});

// Project.rename_group of MDATK 0.1.1 (config.py): renames group `oldName` to `newName` (stripped)
// and every reference in the project, planner steps included: a field holding a group name
// (analysis groups, representation groups, groups shown) and a selection that is the name, or
// comma-separated group names, change when they are `oldName`; analysis instances in data keys
// ("contacts@Aptamer+Target/atom_pairs") and group:<old> tokens in any other text change too, while
// group:<old>er, group:<old>_2 and the notes stay. Returns the changed fields as paths
// ("analyses.contacts.group_a", "layout.panels.1.series.0"; the group's own name not counted);
// nothing for the same name. Unknown keys are not touched (written back unchanged, D9). Throws
// std::out_of_range for an unknown group, ValueError for an invalid or taken name or a name
// several groups have.
std::vector<std::string> renameGroup(Project& project, std::string_view oldName, std::string_view newName);

// Project.group_references of MDATK 0.1.1: the fields that refer to group `name` by the rules
// of renameGroup.
std::vector<std::string> groupReferences(const Project& project, std::string_view name);

// ---- paths -----------------------------------------------------------------------------------

// Rewrites the relative file paths of the project (the paths resolvePaths makes absolute, in its
// order) from `fromFolder` to `toFolder` so that they name the same files: "inputs/a.parm7" in
// /p becomes "../p/inputs/a.parm7" in /q. Empty, absolute and "~" paths stay. Where no relative
// path exists (Windows: another drive or share) the path becomes absolute. Planner files keep their
// ';'-separated layout (each piece rebased); planner settings named by the step type's
// pathSettings are rebased when they are non-empty strings. Returns the changed fields.
std::vector<std::string> rebasePaths(Project& project, std::string_view fromFolder, std::string_view toFolder);

// os.path.relpath(path, start) for absolute, normalised paths of this platform's kind (posixpath;
// on Windows generic ntpath form, nullopt across drives or shares). "." for the same folder.
std::optional<std::string> relativePath(std::string_view path, std::string_view start);

// ---- recovery copies (owner decision D6) -----------------------------------------------------

// A recovery copy is a project file (it opens as a project) with two leading entries,
// recovery_of (the project file, "" for a new project) and recovery_folder (the folder its relative
// paths are relative to), written atomically; the project file itself is never written.

// <mdbridgeConfigDir()>/recovery
std::filesystem::path recoveryDir();

// The recovery file of a document in `dir`: recovery-<16 hex digits of the hash of the absolute
// project path>.mdatk.toml, or recovery-untitled-<process id>-<n>.mdatk.toml for a new project
// (n: the document's number within the process).
std::filesystem::path recoveryPath(const ProjectDocument& document, const std::filesystem::path& dir = recoveryDir());

// Writes the document's recovery copy (creating `dir`); returns its path. Throws
// std::filesystem::filesystem_error.
std::filesystem::path writeRecovery(const ProjectDocument& document, const std::filesystem::path& dir = recoveryDir());

struct RecoveryFile {
    std::filesystem::path file;      // the recovery copy
    std::filesystem::path original;  // the project file it belongs to ("" for a new project)
    std::string folder;              // where its relative paths point from
    DateTime saved;                  // when it was written
};

// The recovery copies in `dir`, newest first (unreadable files are skipped). Never throws.
std::vector<RecoveryFile> findRecoveries(const std::filesystem::path& dir = recoveryDir());

// Opens a recovery copy into `document`: the project of the copy, path() its original (paths
// rebased from the copy's folder when that differs from the original's), modified. Throws
// ProjectError.
void openRecovery(ProjectDocument& document, const RecoveryFile& recovery,
                  std::vector<std::string>* warnings = nullptr);

// Deletes a recovery copy (once saved or discarded). Never throws; false when it could not.
bool removeRecovery(const std::filesystem::path& file);

}  // namespace mdbridge::core
