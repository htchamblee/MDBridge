#include "gui/pages/page_registry.hpp"

#include <functional>
#include <utility>

#include <QMetaObject>
#include <QWidget>

#include "gui/forms/form_context.hpp"
#include "gui/pages/outputs_page.hpp"
#include "gui/pages/placeholder_page.hpp"
#include "gui/pages/structure_page.hpp"

namespace mdbridge::app::pages {

namespace {

struct PageEntry {
    const char* id;
    std::function<QWidget*(forms::FormContext&)> make;
};

QWidget* placeholder(forms::FormContext& c, const char* id) { return new PlaceholderPage(c, QString::fromLatin1(id)); }

// MDATK's order (app.py:226-254). A page not delivered yet is a PlaceholderPage (every field of
// the page from the model) until its package replaces the line: input, groups, analyses (W1-B),
// layout (W2-A).
const std::vector<PageEntry>& entries() {
    static const std::vector<PageEntry> list{
        {"input", [](forms::FormContext& c) { return placeholder(c, "input"); }},
        {"groups", [](forms::FormContext& c) { return placeholder(c, "groups"); }},
        {"analyses", [](forms::FormContext& c) { return placeholder(c, "analyses"); }},
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
