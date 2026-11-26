#pragma once
// The settings pages of the analysis window, in MDATK's order (core::kProblemPages: input, groups,
// analyses, structure, layout, outputs). Each page package adds its page to createPages(); the
// window shows them in its side bar and goes to a problem's field with revealField.
//
// A page is a QWidget with pageId() (one of core::kProblemPages' ids), title() and
// revealField(path) (scrolls to the field, opens its section, focuses it; false when the field is
// not on the page). The helpers below call them through Qt's meta-object system, so a page need not
// derive from FormPage as long as the three are Q_INVOKABLE.

#include <vector>

#include <QString>

class QWidget;

namespace mdbridge::app::forms {
class FormContext;
}

namespace mdbridge::app::pages {

// The pages there are, in MDATK's order, bound to `context` (which must outlive them). The caller
// owns them (they have no parent).
std::vector<QWidget*> createPages(forms::FormContext& context);

// The ids of the pages createPages() makes, in order.
std::vector<QString> pageIds();

QString pageIdOf(QWidget* page);
QString pageTitleOf(QWidget* page);
bool revealFieldOn(QWidget* page, const QString& path);
// The page is a PlaceholderPage (its own page is not there yet; placeholder_page.hpp).
bool isPlaceholder(QWidget* page);

}  // namespace mdbridge::app::pages
