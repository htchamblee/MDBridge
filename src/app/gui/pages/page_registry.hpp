#pragma once
// The settings pages of the side bar (milestone 5 wave 1), in MDATK's order (core::kProblemPages):
// Input & cleaning, Residue groups, Analyses, Structure view, Layout & panels, Outputs. Each page is
// a QWidget with pageId(), title() and revealField(path) (invokable, see settings_page.hpp); pages
// of later packages are added here as they land.

#include <vector>

class QWidget;

namespace mdbridge::app::forms {
class FormContext;
}

namespace mdbridge::app::pages {

// New pages bound to `context`, in MDATK's order; the caller owns them (no parent set).
std::vector<QWidget*> createPages(forms::FormContext& context);

}  // namespace mdbridge::app::pages
