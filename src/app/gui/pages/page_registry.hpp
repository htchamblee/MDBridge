#pragma once
// The settings pages of the main window's side bar, in MDATK's order (core::kProblemPages: Input &
// cleaning, Residue groups, Analyses, Structure view, Layout & panels, Outputs).
//
// A page is a QWidget class XPage(forms::FormContext&, QWidget* parent) with QString pageId(),
// QString title() and bool revealField(const QString& path) (scrolls to the field and focuses it;
// false when the field is not on the page), and sections that are PageSection
// (page_section.hpp, objectName "section:<page>.<section>"). The window puts each page into a
// scroll area of the side bar's fixed width (pages scroll vertically, never sideways), so a page
// has no scroll area of its own.
//
// createPages() is the one place that lists them: each package adds its page by replacing the
// placeholder's line (`entry(new InputPage(context))`).

#include <functional>
#include <vector>

#include <QString>

class QWidget;

namespace mdbridge::app::forms {
class FormContext;
}

namespace mdbridge::app::pages {

struct PageEntry {
    QString id;
    QString title;
    QWidget* widget = nullptr;  // the caller takes ownership (parents it)
    std::function<bool(const QString& path)> revealField;
    bool placeholder = false;   // a PlaceholderPage (the page's own is not there yet)
};

// Any page class of the form above.
template <class Page>
PageEntry entry(Page* page, bool placeholder = false) {
    return PageEntry{page->pageId(), page->title(), page,
                     [page](const QString& path) { return page->revealField(path); }, placeholder};
}

// The six pages, in MDATK's order, bound to `context`.
std::vector<PageEntry> createPages(forms::FormContext& context);

}  // namespace mdbridge::app::pages
