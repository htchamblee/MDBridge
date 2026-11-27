#include "gui/tools/tool_windows.hpp"

#include <utility>

#include "gui/main_window.hpp"
#include "gui/tools/builder_window.hpp"
#include "gui/tools/pdb_tools_window.hpp"

namespace mdbridge::app::tools {

// ---- MainWindowHost --------------------------------------------------------------------------

QWidget* MainWindowHost::window() { return &m_window; }

core::ProjectDocument& MainWindowHost::document() { return m_window.document(); }

std::filesystem::path MainWindowHost::projectPath() const { return m_window.projectFile(); }

void MainWindowHost::appendLog(const QString& text) { m_window.appendLog(text); }

void MainWindowHost::error(const QString& text, const QString& detail) { m_window.showError(text, detail); }

void MainWindowHost::info(const QString& text, const QString& detail) { m_window.showInfo(text, detail); }

bool MainWindowHost::saveProject(bool saveAs) { return saveAs ? m_window.saveAs() : m_window.save(); }

Prompter& MainWindowHost::prompter() { return m_window.prompter(); }

// ---- ToolWindows -----------------------------------------------------------------------------

ToolWindows::ToolWindows(std::unique_ptr<Host> host) : m_host(std::move(host)) {}

ToolWindows::~ToolWindows() {
    delete m_pdbTools.data();
    delete m_builder.data();  // a running build is cancelled and waited for (its JobRunner)
}

namespace {

template <class W>
W* present(QPointer<W>& slot, Host& host) {
    if (!slot) {
        slot = new W(host, host.window());
    }
    slot->show();
    slot->raise();
    slot->activateWindow();
    return slot;
}

}  // namespace

PdbToolsWindow* ToolWindows::openPdbTools() { return present(m_pdbTools, *m_host); }

BuilderWindow* ToolWindows::openBuilder() { return present(m_builder, *m_host); }

PdbToolsWindow* ToolWindows::pdbTools() const { return m_pdbTools.data(); }

BuilderWindow* ToolWindows::builder() const { return m_builder.data(); }

}  // namespace mdbridge::app::tools
