#include "gui/pages/page_registry.hpp"

#include "gui/forms/form_context.hpp"
#include "gui/pages/placeholder_page.hpp"

namespace mdbridge::app::pages {

namespace {

PageEntry placeholder(forms::FormContext& context, const char* id) {
    return entry(new PlaceholderPage(context, QString::fromLatin1(id)), true);
}

}  // namespace

std::vector<PageEntry> createPages(forms::FormContext& context) {
    std::vector<PageEntry> out;
    out.push_back(placeholder(context, "input"));      // W1-B: InputPage
    out.push_back(placeholder(context, "groups"));     // W1-B: GroupsPage
    out.push_back(placeholder(context, "analyses"));   // W1-B: AnalysesPage
    out.push_back(placeholder(context, "structure"));  // W1-C: StructurePage
    out.push_back(placeholder(context, "layout"));     // W2-A: LayoutPage
    out.push_back(placeholder(context, "outputs"));    // W1-C: OutputsPage
    return out;
}

}  // namespace mdbridge::app::pages
