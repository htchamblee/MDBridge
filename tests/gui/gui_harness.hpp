#pragma once
// The GUI test harness (tests/data/README.md, "gui"): action scripts run against MainWindow, the
// files written compared with tests/data/gui, widgets found by objectName (field widgets use their
// project key), screenshots for CI artifacts (never compared).

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <QString>

#include "gui/dialogs/prompts.hpp"
#include "gui/forms/field_widgets.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/value.hpp"

class QWidget;

namespace mdbridge::app {
class MainWindow;
}

namespace mdbridge::guitest {

// A step of a case failed or a case file is malformed; the message says which and why.
class CaseError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct GuiCase {
    std::filesystem::path source;  // the case file
    std::string name;              // its path under the cases folder, without .json ("harness/open_save_default")
    std::string description;
    std::optional<core::DateTime> clock;
    std::vector<std::pair<std::string, std::string>> files;     // case-folder name, path under the data folder
    core::Array actions;                                         // tables with "action"
    std::vector<std::pair<std::string, std::string>> expected;  // case-folder name, path under the data folder
    std::vector<std::string> needs;                              // "needs": features beyond this harness (W0-A cases)
    std::string start;                       // "start.project": the file the window starts with ("" none)
    std::string settings;                    // "settings": settings.toml before the start (path under tests/data)
    std::optional<std::vector<std::string>> log;  // "log": the log lines expected
    core::Value outcome;                     // "outcome" (a table; Null without one)
    std::optional<std::pair<int, int>> canvas;  // "canvas": the figure canvas's size in pixels

    bool needsAny(std::string_view need) const;
};

// The entries of `c.needs` this harness does not support yet (it runs "type" and "table"; a case
// that needs anything else is skipped, tests/data/README.md "gui").
std::vector<std::string> unsupportedNeeds(const GuiCase& c);

// Every *.json under `dir` (recursively), sorted.
std::vector<std::filesystem::path> findCases(const std::filesystem::path& dir);

// Reads a case file; `name` is its path relative to `casesDir`. Throws CaseError.
GuiCase loadCase(const std::filesystem::path& file, const std::filesystem::path& casesDir);

// The project file clock fixed at `start`, one second further on every reading, while it lives.
class FixedClock {
public:
    explicit FixedClock(const core::DateTime& start);
    ~FixedClock();
    FixedClock(const FixedClock&) = delete;
    FixedClock& operator=(const FixedClock&) = delete;

private:
    core::Clock m_previous;
};

// The one widget under `root` (itself included) whose objectName is `key`. Throws CaseError when
// there is none or more than one.
QWidget* findField(QWidget* root, const QString& key);

// Enters `value` into the widget as a user would (tests/data/README.md, action "set"): drop-downs
// choose the item whose project value (item data) is `value` (W0-A's cases), editable ones and text
// fields get it typed (select all, type, Enter; a list of strings typed as "a, b"), check boxes are
// toggled with Space, sliders moved to the value (times their "mdbridgeScale" property). A list
// value entered into a widget that is none of these goes to its parts: check boxes with an
// "mdbridgeChoice" property are ticked for the values listed and unticked for the others (unticks
// first, then ticks in the list's order: output.charts.formats), else item k goes to the child
// named "<key>.k" (the rotation sliders "visual.rotation.0" ...). Throws CaseError for widget types
// and values it cannot enter.
void setField(QWidget* widget, const core::Value& value);

// Replaces the widget's text with `text` and commits it (action "type"): select all, type (""
// deletes the selection, clearing the field), Enter; multi-line text is committed as on leaving
// the field. Throws CaseError.
void typeField(QWidget* widget, const QString& text);

// "set" / "type" by key under `root`: the widget named `key`, else the list editor cell of that
// field path ("groups.3.name", forms::findCell): the cell's editor opened (QTableView::edit), the
// value entered there as above and Enter pressed; check cells are toggled with Space. Throws
// CaseError.
void setKey(QWidget* root, const QString& key, const core::Value& value);
void typeKey(QWidget* root, const QString& key, const QString& text);

// Action "table": in the list editor of `key` (its field path), "add" clicks <key>:add, "remove"
// selects row `row` and clicks <key>:remove, "move" selects row `row` and clicks <key>:up (by < 0)
// or <key>:down (by > 0), |by| times. Throws CaseError.
void tableAction(QWidget* root, const QString& key, const std::string& what, int row, int by);

// Saves `widget` as <dir>/<name>.png (the folder is made); returns the file. Throws CaseError.
std::filesystem::path saveScreenshot(QWidget* widget, const std::filesystem::path& dir, const std::string& name);

// The application font of the GUI tests: the bundled sans face (render::FontFace::Sans) at Qt's
// default 9 pt, so text widths are the same on every machine (Qt's offscreen platform on Windows
// finds no font at all and draws boxes; Linux CI and developers' machines have different fonts)
// and layout tests such as WindowTests::pagesFitTheSideBar measure what CI measures. Call after
// the QApplication exists; a warning (and Qt's default font) when the bundled fonts cannot be
// registered.
void useBundledFont();

// The screenshot folder: $MDBRIDGE_GUI_SCREENSHOTS, else `fallback`.
std::filesystem::path screenshotDir(const std::filesystem::path& fallback);

// `text` with every occurrence of the folder `folder` replaced by `name`. With `anySeparator` (on
// Windows) '\\' and '/' match each other in the folder, and the separators of the path that follows
// it, up to a space, quote, colon, comma, semicolon or closing bracket, become '/'.
std::string replaceFolder(std::string text, const std::string& folder, const std::string& name, bool anySeparator);

// The text as compared: CR LF -> LF, the case folder's absolute path (generic and native forms;
// on Windows with either separator, the path after it with '/') -> {case_dir}.
std::string normalise(std::string text, const std::filesystem::path& caseDir);

// "" when equal, else where the texts first differ ("line 12: expected ... got ...").
std::string compareText(const std::string& actual, const std::string& expected);

// The window's questions, messages and terminal lines during a case: each recorded with the index
// of the action that caused it (-1: the start), questions answered from the action's "answers"
// (a question without an answer gets the first choice, tests/data/README.md "gui").
class RecordingPrompter : public app::Prompter {
public:
    QString ask(QWidget* parent, const app::Question& question) override;
    void message(QWidget* parent, app::MessageKind kind, const QString& text, const QString& detail) override;
    void printError(const QString& line) override;

    int action = -1;
    core::Value answers;                  // the current action's "answers" (a table, or Null)
    core::Array questions;                // {"action", "kind", "path", "answer"}
    core::Array dialogs;                  // {"action", "kind", "text", "detail"}
    std::vector<std::string> stderrLines;
};

// Runs a case's actions against `window`, files in `caseDir`, screenshots to `screenshots`
// (named <case name with / as _>-<name>.png). Throws CaseError naming the action.
class CaseRunner {
public:
    CaseRunner(app::MainWindow& window, std::filesystem::path caseDir, std::filesystem::path screenshots,
               std::string caseName);

    ~CaseRunner();
    CaseRunner(const CaseRunner&) = delete;
    CaseRunner& operator=(const CaseRunner&) = delete;

    // The window started with `project` (case folder; as `mdbridge gui PROJECT`).
    void start(const std::string& project);
    void run(const core::Array& actions);
    void runAction(const core::Value& action);
    // "open" may fail (the case needs "open_failure"): its error is a recorded dialog then.
    void setOpenMayFail(bool may) { m_openMayFail = may; }
    // MDATK's canvas records of the case (outcome.canvas: hit rectangles and the pointer positions
    // each canvas action resolved to): the canvas actions send their events at those positions, so
    // that the drags are MDATK's mouse paths whatever this canvas's geometry; without a record the
    // points are resolved from this canvas's hit areas.
    void setCanvasRecords(const core::Value& records) { m_canvasRecords = records; }
    // Where this canvas's hit areas differ from MDATK's by more than a few pixels (one line each).
    const std::vector<std::string>& geometryNotes() const { return m_geometryNotes; }

    const std::vector<std::filesystem::path>& screenshots() const { return m_shots; }
    const RecordingPrompter& prompter() const { return m_prompter; }
    const core::Array& snapshots() const { return m_snapshots; }

private:
    void waitUntilLoaded(int timeoutMs);
    void waitForPreview();
    void snapshot(const core::Value& action);
    void canvasAction(const core::Value& action, const std::string& what);
    QStringList dialogAnswer(const app::forms::FileDialogRequest& request) const;

    app::MainWindow& m_window;
    RecordingPrompter m_prompter;
    core::Array m_snapshots;
    core::Value m_current;  // the action being run
    int m_index = -1;
    bool m_openMayFail = false;
    core::Value m_canvasRecords;
    std::vector<std::string> m_geometryNotes;
    app::forms::FileDialogFunction m_previousDialog;
    std::filesystem::path m_caseDir;
    std::filesystem::path m_screenshotDir;
    std::string m_caseName;
    std::vector<std::filesystem::path> m_shots;
};

// Copies the case's files into `caseDir` (which must exist). Throws CaseError.
void prepareCase(const GuiCase& c, const std::filesystem::path& dataDir, const std::filesystem::path& caseDir);

// Every mismatch of the case's expected files, one per line; "" when all are equal.
std::string checkExpected(const GuiCase& c, const std::filesystem::path& dataDir, const std::filesystem::path& caseDir);

// Every mismatch of the case's outcome (tests/data/README.md "gui": log, questions, dialogs,
// snapshots, stderr, settings) with what `runner` recorded, `log` (the window's log lines) and the
// settings file; texts with {case_dir} and {config_dir} for the two folders. "" when all agree.
std::string checkOutcome(const GuiCase& c, const CaseRunner& runner, const std::vector<std::string>& log,
                         const std::filesystem::path& dataDir, const std::filesystem::path& caseDir,
                         const std::filesystem::path& configDir);

// The text with the case and config folders (generic and native forms; on Windows with either
// separator, the paths after them with '/') as {case_dir} / {config_dir}.
std::string placeholders(std::string text, const std::filesystem::path& caseDir, const std::filesystem::path& configDir);

// The words of the TOML parser in "not valid TOML: <words> (at line ..." replaced by a marker.
std::string maskParserWords(std::string text);

// Two values as the snapshots compare them: numbers within 1e-9 relative (NaN equals NaN; an
// integer equals the same float), tables by key whatever their order, everything else exactly.
bool sameValue(const core::Value& a, const core::Value& b);

}  // namespace mdbridge::guitest
