#pragma once
// The workflow tools of a host window (milestone 6 W3-B, owner decision D1): the PDB clean-up and
// the structure builder, opened from the analysis window's Workflow menu as in MDATK v0.0.4.1
// (app.py open_pdb_tools / open_builder: one window of each, shown again when it is open), each
// knowing its host only through tools::Host (MDATK 0.1.1's gui/host.py). Milestone 7 moves them
// into the `rnado` program by giving them its host; the planner and the script wizard join them
// then.
//
// MainWindowHost is the analysis window as such a host.

#include <filesystem>
#include <memory>

#include <QPointer>

#include "gui/tools/host.hpp"

namespace mdbridge::app {
class MainWindow;
}

namespace mdbridge::app::tools {

class BuilderWindow;
class PdbToolsWindow;

// The analysis window as the tools' host: its project file (none for an unsaved project or a
// parameter log), log, messages, save and questions.
class MainWindowHost : public Host {
public:
    explicit MainWindowHost(MainWindow& window) : m_window(window) {}
    QWidget* window() override;
    core::ProjectDocument& document() override;
    std::filesystem::path projectPath() const override;
    void appendLog(const QString& text) override;
    void error(const QString& text, const QString& detail = QString()) override;
    void info(const QString& text, const QString& detail = QString()) override;
    bool saveProject(bool saveAs = false) override;
    Prompter& prompter() override;

private:
    MainWindow& m_window;
};

class ToolWindows {
public:
    // The windows are children of host->window(); they are closed (deleted) with this object,
    // before the host.
    explicit ToolWindows(std::unique_ptr<Host> host);
    ~ToolWindows();
    ToolWindows(const ToolWindows&) = delete;
    ToolWindows& operator=(const ToolWindows&) = delete;

    Host& host() const { return *m_host; }

    // Opens the window (a new one, or the open one shown and raised) and returns it.
    PdbToolsWindow* openPdbTools();
    BuilderWindow* openBuilder();
    // The open windows (null once closed).
    PdbToolsWindow* pdbTools() const;
    BuilderWindow* builder() const;

private:
    std::unique_ptr<Host> m_host;
    QPointer<PdbToolsWindow> m_pdbTools;
    QPointer<BuilderWindow> m_builder;
};

}  // namespace mdbridge::app::tools
