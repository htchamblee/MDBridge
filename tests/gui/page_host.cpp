#include "page_host.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>

#include <QMetaObject>
#include <QStackedWidget>
#include <QVBoxLayout>

#include "gui/forms/form_context.hpp"
#include "gui/pages/page_registry.hpp"
#include "mdbridge/core/paths.hpp"

namespace mdbridge::guitest {

PageHost::PageHost() {
    m_context = std::make_unique<app::forms::FormContext>(m_document);
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_stack = new QStackedWidget(this);
    layout->addWidget(m_stack);
    m_pages = app::pages::createPages(*m_context);
    for (QWidget* page : m_pages) {
        m_stack->addWidget(page);
    }
    resize(520, 900);  // MDATK's side bar is 520 px wide
}

PageHost::~PageHost() {
    // the pages' field widgets go before the context they are registered with
    for (QWidget* page : m_pages) {
        delete page;
    }
    m_pages.clear();
}

bool PageHost::openProject(const std::filesystem::path& path) {
    try {
        m_document.open(path);
    } catch (const std::exception& e) {
        std::cerr << "cannot open " << core::pathToUtf8(path) << ": " << e.what() << "\n";
        return false;
    }
    return true;
}

bool PageHost::saveProject(const std::filesystem::path& path) {
    try {
        if (path.empty()) {
            m_document.save();
        } else {
            m_document.saveAs(path);
        }
    } catch (const std::exception& e) {
        std::cerr << "cannot save: " << e.what() << "\n";
        return false;
    }
    return true;
}

void PageHost::commitPendingEdits() { m_context->commitPendingEdits(); }

void PageHost::showPage(QWidget* page) { m_stack->setCurrentWidget(page); }

void PageHost::reveal(const QString& key) {
    for (QWidget* page : m_pages) {
        bool found = false;
        m_stack->setCurrentWidget(page);  // revealField scrolls within the page shown
        QMetaObject::invokeMethod(page, "revealField", Q_RETURN_ARG(bool, found), Q_ARG(QString, key));
        if (found) {
            return;
        }
    }
}

QString PageHost::currentPage() const {
    QString id;
    if (QWidget* page = m_stack->currentWidget()) {
        QMetaObject::invokeMethod(page, "pageId", Q_RETURN_ARG(QString, id));
    }
    return id;
}

bool runsInPageHost(const GuiCase& c) {
    // until the window (W1-A) shows the pages: the cases of the pages this host has
    static const std::vector<std::string> cases{"pages/input_page", "pages/groups_page", "pages/analyses_page"};
    return std::find(cases.begin(), cases.end(), c.name) != cases.end();
}

const std::vector<std::string>& pageHostNeeds() {
    static const std::vector<std::string> needs{"lossless", "snapshot", "start"};
    return needs;
}

}  // namespace mdbridge::guitest
