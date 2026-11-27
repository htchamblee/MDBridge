#pragma once
// What the workflow tools need from the window that hosts them (milestone 6 W3-B, owner decision
// D1): MDATK 0.1.1's gui/host.py `Host`. The PDB clean-up and the structure builder (and, with
// milestone 7, the planner and the script wizard) use nothing else of the analysis window, so that
// milestone 7 can move them into the separate `rnado` program by giving them another host; a test
// can host them with any object that implements this.
//
// MDATK 0.1.1's host also runs and stops the planner's scripts (run_scripts, stop_scripts,
// stop_step), hands a plan's files to the analysis window (use_for_analysis) and adjusts the plan
// without an edit (without_edit): milestone 7 adds those with the planner.

#include <filesystem>

#include <QString>

class QWidget;

namespace mdbridge::core {
class ProjectDocument;
}

namespace mdbridge::app {

class Prompter;

namespace tools {

class Host {
public:
    Host() = default;
    virtual ~Host() = default;
    Host(const Host&) = delete;
    Host& operator=(const Host&) = delete;

    // The window the tools stay above (Qt's transient parent); null for none.
    virtual QWidget* window() = 0;
    // The open project (0.1.1: `project`; the planner and the wizard of milestone 7 edit it).
    virtual core::ProjectDocument& document() = 0;
    // The project file whose step log records what the tools run (0.1.1: `project_path`); empty
    // for an unsaved project or an opened parameter log, whose steps are not logged.
    virtual std::filesystem::path projectPath() const = 0;
    // A message for the window's log ("PDB clean-up: wrote out.pdb").
    virtual void appendLog(const QString& text) = 0;
    // A problem / a result: logged and shown.
    virtual void error(const QString& text, const QString& detail = QString()) = 0;
    virtual void info(const QString& text, const QString& detail = QString()) = 0;
    // Saves the project (asks for a file without one, or always with `saveAs`); false when it was
    // not saved.
    virtual bool saveProject(bool saveAs = false) = 0;
    // The questions and messages of the host (tests answer them).
    virtual Prompter& prompter() = 0;
};

}  // namespace tools
}  // namespace mdbridge::app
