#pragma once
// A settings page built from the project model until its own page is delivered (milestone 5 W1-A;
// the pages of W1-B and W2-A replace it in page_registry.cpp): every field the page shows
// (core::problemPage) with the widget forms::createFieldWidget gives it, one folding section per
// part of the project ("inputs", "analyses.contacts", "output.video" ...), lists of records as list
// editors. It edits the document as the real page will, so the window is usable and the GUI cases
// that only need fields run, but it has none of MDATK's page layout, hints or buttons.

#include <QString>

#include "gui/pages/form_page.hpp"

namespace mdbridge::app::pages {

class PlaceholderPage : public FormPage {
    Q_OBJECT

public:
    // `id`: one of core::kProblemPages.
    PlaceholderPage(forms::FormContext& context, const QString& id, QWidget* parent = nullptr);
};

}  // namespace mdbridge::app::pages
