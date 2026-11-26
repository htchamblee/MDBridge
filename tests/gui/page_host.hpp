#pragma once
// A test-only host of the settings pages (milestone 5 W1-B): a document, its FormContext and the
// pages of pages::createPages() in a stacked widget, the page of a field shown before each action
// on it, as the window's side bar will. It lets the per-page GUI cases (tests/data/gui/pages) run
// before the window (W1-A) has the pages; once it has, the cases run against the window instead.

#include <filesystem>
#include <memory>
#include <vector>

#include <QWidget>

#include "gui_harness.hpp"
#include "mdbridge/core/project_document.hpp"

class QStackedWidget;

namespace mdbridge::app::forms {
class FormContext;
}

namespace mdbridge::guitest {

class PageHost : public QWidget, public CaseHost {
public:
    PageHost();
    ~PageHost() override;

    QWidget* widget() override { return this; }
    bool openProject(const std::filesystem::path& path) override;
    bool saveProject(const std::filesystem::path& path) override;
    bool settled() const override { return true; }
    core::ProjectDocument* document() override { return &m_document; }
    void commitPendingEdits() override;
    void reveal(const QString& key) override;
    QString currentPage() const override;

    const std::vector<QWidget*>& pages() const { return m_pages; }
    void showPage(QWidget* page);
    app::forms::FormContext& context() { return *m_context; }

private:
    core::ProjectDocument m_document;
    std::unique_ptr<app::forms::FormContext> m_context;
    QStackedWidget* m_stack = nullptr;
    std::vector<QWidget*> m_pages;
};

// The cases this host runs (the pages it has), and the needs it supports besides the harness's.
bool runsInPageHost(const GuiCase& c);
const std::vector<std::string>& pageHostNeeds();

}  // namespace mdbridge::guitest
