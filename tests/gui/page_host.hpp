#pragma once
// A test-only host for the settings pages (milestone 5 W1-C): until the analysis window has its
// pages (W1-A), the GUI cases whose actions are all on pages that exist run against this host
// instead of MainWindow. It holds a ProjectDocument (the case's `start` project opened, as
// `mdbridge gui PATH` does), its FormContext and every page of pages::createPages() in a tab
// widget, and runs the page-level actions of tests/data/README.md "gui": set, type, table,
// click, toggle, open, save, snapshot (keys and density_status) and screenshot. The snapshots are
// compared with the case's outcome.snapshots.
//
// Not covered (window-level; the case runs here only when it needs none of them): the recent-
// projects file (outcome.settings), the log lines, dialogs and questions, canvas actions, New /
// Close / Save as, the page shown, the modified marker and slot editors.

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <QWidget>

#include "gui_harness.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/value.hpp"

class QTabWidget;

namespace mdbridge::app::forms {
class FormContext;
}

namespace mdbridge::guitest {

class PageHost : public QWidget {
    Q_OBJECT

public:
    PageHost();
    ~PageHost() override;

    core::ProjectDocument& document() { return m_document; }
    app::forms::FormContext& context() { return *m_context; }
    const std::vector<QWidget*>& pages() const { return m_pages; }
    QWidget* page(const QString& id) const;
    QTabWidget* tabs() const { return m_tabs; }

private:
    core::ProjectDocument m_document;
    std::unique_ptr<app::forms::FormContext> m_context;
    QTabWidget* m_tabs = nullptr;
    std::vector<QWidget*> m_pages;
};

// Whether the host can run case `c`: a start project, needs only of {start, lossless, type, table,
// snapshot, toggle}, and every key an action names is on a page the host has. `why` gets the
// reason when not.
bool pageHostRuns(const GuiCase& c, std::string& why);

// Runs `c` against a host (files prepared in `caseDir`, screenshots of each page shown by the
// case to `screenshots`); returns every mismatch ("" when the case passes), the expected files
// included.
std::string runOnPageHost(const GuiCase& c, const std::filesystem::path& dataDir, const std::filesystem::path& caseDir,
                          const std::filesystem::path& screenshots, std::vector<std::filesystem::path>* shots = nullptr);

// "" when `actual` equals `expected` (tables by key, numbers within 1e-9 relative, an integer
// equal to a float of the same value), else the first difference, named by `where`.
std::string compareValues(const core::Value& actual, const core::Value& expected, const std::string& where);

}  // namespace mdbridge::guitest
