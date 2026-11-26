#include "gui/pages/placeholder_page.hpp"

#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "gui/forms/field_factory.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

namespace {

using forms::toQString;

// The part of the project a field belongs to: "analyses.contacts" and "output.video" for the
// fields of an analysis or an output, else the first key ("inputs", "groups", "notes").
std::string sectionPath(const std::string& path) {
    const std::size_t first = path.find('.');
    if (first == std::string::npos) return path;
    const std::string head = path.substr(0, first);
    if (head == "analyses" || head == "output") {
        const std::size_t second = path.find('.', first + 1);
        if (second != std::string::npos) return path.substr(0, second);
    }
    return head;
}

// The last part of a described location ("Analyses › Pairwise contacts" -> "Pairwise contacts").
QString lastPart(const std::string& location) {
    const QString text = toQString(location);
    const qsizetype at = text.lastIndexOf(QStringLiteral(" › "));
    return at < 0 ? text : text.mid(at + 3);
}

bool inside(const std::string& path, const std::string& list) {
    return path.size() > list.size() && path.compare(0, list.size(), list) == 0 && path[list.size()] == '.';
}

// The fields of `page`: list editors for lists of records (their items are its cells).
std::vector<std::string> pagePaths(const core::Project& project, const std::string& page) {
    std::vector<std::string> out;
    std::vector<std::string> lists;
    for (const std::string& path : forms::formFieldPaths(project)) {
        if (core::problemPage(path) != page) continue;
        if (std::any_of(lists.begin(), lists.end(), [&](const std::string& list) { return inside(path, list); })) {
            continue;
        }
        const std::optional<core::FieldInfo> info = core::fieldInfo(project, path);
        if (!info || info->kind == core::FieldKind::Record) continue;
        if (info->kind == core::FieldKind::Records || info->kind == core::FieldKind::NumberRows) {
            lists.push_back(path);
        }
        out.push_back(path);
    }
    return out;
}

}  // namespace

PlaceholderPage::PlaceholderPage(forms::FormContext& context, const QString& id, QWidget* parent)
    : FormPage(context, id, toQString(core::problemPageTitle(id.toStdString())), parent) {
    const core::Project& project = context.document().project();
    std::map<std::string, Section*> sections;
    for (const std::string& path : pagePaths(project, id.toStdString())) {
        const std::string part = sectionPath(path);
        Section*& section = sections[part];
        if (!section) {
            const std::size_t dot = part.rfind('.');
            QString sectionId = toQString(dot == std::string::npos ? part : part.substr(dot + 1));
            QString heading = lastPart(core::describeField(part, &project));
            if (part == "analyses" || part == "output") {  // the fields directly under analyses / output
                sectionId = QStringLiteral("general");
                heading = part == "analyses" ? tr("General") : tr("Locations");
            }
            section = addSection(sectionId, heading);
            if (sections.size() == 1) {
                section->form()->addHint(tr("Every setting of this page, one field each; the page's own layout "
                                            "comes with a later version."));
            }
        }
        forms::FieldWidget* widget = forms::createFieldWidget(context, path, section->form());
        if (!widget) continue;
        QString label = widget->label();
        if (label.isEmpty()) label = lastPart(core::describeField(path, &project));
        section->form()->add(label, widget);
    }
}

}  // namespace mdbridge::app::pages
