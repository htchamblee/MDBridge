#include "gui/pages/page_registry.hpp"

#include "gui/pages/analyses_page.hpp"
#include "gui/pages/groups_page.hpp"
#include "gui/pages/input_page.hpp"

namespace mdbridge::app::pages {

std::vector<QWidget*> createPages(forms::FormContext& context) {
    return {
        new InputPage(context),
        new GroupsPage(context),
        new AnalysesPage(context),
        // Structure view, Layout & panels, Outputs: W1-C, W2-A
    };
}

}  // namespace mdbridge::app::pages
