#pragma once
// The GUI test harness (tests/data/README.md, "gui"): action scripts run against MainWindow, the
// files written compared with tests/data/gui, widgets found by objectName (field widgets use their
// project key), screenshots for CI artifacts (never compared).

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <QString>

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
};

// The entries of `c.needs` this harness does not support yet (it supports none of them: a case that
// needs any is skipped, tests/data/README.md "gui").
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

// Enters `value` into the widget as a user would (tests/data/README.md, action "set"). Throws
// CaseError for widget types and values it cannot enter.
void setField(QWidget* widget, const core::Value& value);

// Saves `widget` as <dir>/<name>.png (the folder is made); returns the file. Throws CaseError.
std::filesystem::path saveScreenshot(QWidget* widget, const std::filesystem::path& dir, const std::string& name);

// The screenshot folder: $MDBRIDGE_GUI_SCREENSHOTS, else `fallback`.
std::filesystem::path screenshotDir(const std::filesystem::path& fallback);

// The text as compared: CR LF -> LF, the case folder's absolute path (generic form) -> {case_dir}.
std::string normalise(std::string text, const std::filesystem::path& caseDir);

// "" when equal, else where the texts first differ ("line 12: expected ... got ...").
std::string compareText(const std::string& actual, const std::string& expected);

// Runs a case's actions against `window`, files in `caseDir`, screenshots to `screenshots`
// (named <case name with / as _>-<name>.png). Throws CaseError naming the action.
class CaseRunner {
public:
    CaseRunner(app::MainWindow& window, std::filesystem::path caseDir, std::filesystem::path screenshots,
               std::string caseName);

    void run(const core::Array& actions);
    void runAction(const core::Value& action);

    const std::vector<std::filesystem::path>& screenshots() const { return m_shots; }

private:
    void waitUntilLoaded(int timeoutMs);

    app::MainWindow& m_window;
    std::filesystem::path m_caseDir;
    std::filesystem::path m_screenshotDir;
    std::string m_caseName;
    std::vector<std::filesystem::path> m_shots;
};

// Copies the case's files into `caseDir` (which must exist). Throws CaseError.
void prepareCase(const GuiCase& c, const std::filesystem::path& dataDir, const std::filesystem::path& caseDir);

// Every mismatch of the case's expected files, one per line; "" when all are equal.
std::string checkExpected(const GuiCase& c, const std::filesystem::path& dataDir, const std::filesystem::path& caseDir);

}  // namespace mdbridge::guitest
