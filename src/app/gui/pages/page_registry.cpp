#include "gui/pages/page_registry.hpp"

#include <functional>
#include <utility>

#include <QMetaObject>
#include <QWidget>

#include "gui/forms/form_context.hpp"
#include "gui/pages/analyses_page.hpp"
#include "gui/pages/groups_page.hpp"
#include "gui/pages/input_page.hpp"
#include "gui/pages/outputs_page.hpp"
#include "gui/pages/placeholder_page.hpp"
#include "gui/pages/structure_page.hpp"

namespace mdbridge::app::pages {

namespace {

struct PageEntry {
    const char* id;
    std::function<QWidget*(forms::FormContext&)> make;
};

// MDATK's order (app.py:226-254). The Layout & panels page is a PlaceholderPage (every field of
// the page from the model) until W2-A replaces its line.
QWidget* placeholder(forms::FormContext& c, const char* id) { return new PlaceholderPage(c, QString::fromLatin1(id)); }

const std::vector<PageEntry>& entries() {
    static const std::vector<PageEntry> list{
        {"input", [](forms::FormContext& c) -> QWidget* { return new InputPage(c); }},
        {"groups", [](forms::FormContext& c) -> QWidget* { return new GroupsPage(c); }},
        {"analyses", [](forms::FormContext& c) -> QWidget* { return new AnalysesPage(c); }},
        {"structure", [](forms::FormContext& c) -> QWidget* { return new StructurePage(c); }},
        {"layout", [](forms::FormContext& c) { return placeholder(c, "layout"); }},
        {"outputs", [](forms::FormContext& c) -> QWidget* { return new OutputsPage(c); }},
    };
    return list;
}

}  // namespace

std::vector<QWidget*> createPages(forms::FormContext& context) {
    std::vector<QWidget*> out;
    for (const PageEntry& e : entries()) {
        out.push_back(e.make(context));
    }
    return out;
}

bool isPlaceholder(QWidget* page) { return qobject_cast<PlaceholderPage*>(page) != nullptr; }

std::vector<QString> pageIds() {
    std::vector<QString> out;
    for (const PageEntry& e : entries()) {
        out.push_back(QString::fromLatin1(e.id));
    }
    return out;
}

QString pageIdOf(QWidget* page) {
    QString id;
    if (page) {
        QMetaObject::invokeMethod(page, "pageId", Qt::DirectConnection, qReturnArg(id));
    }
    return id;
}

QString pageTitleOf(QWidget* page) {
    QString title;
    if (page) {
        QMetaObject::invokeMethod(page, "title", Qt::DirectConnection, qReturnArg(title));
    }
    return title;
}

bool revealFieldOn(QWidget* page, const QString& path) {
    bool shown = false;
    if (page) {
        QMetaObject::invokeMethod(page, "revealField", Qt::DirectConnection, qReturnArg(shown), path);
    }
    return shown;
}

}  // namespace mdbridge::app::pages
