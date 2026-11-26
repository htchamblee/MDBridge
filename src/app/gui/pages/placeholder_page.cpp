#include "gui/pages/placeholder_page.hpp"

#include <algorithm>
#include <map>
#include <string_view>

#include <QFormLayout>
#include <QLabel>
#include <QVBoxLayout>

#include "gui/forms/field_factory.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/pages/page_section.hpp"
#include "gui/problems/problems_panel.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

namespace {

using forms::toQString;

// The part of the project a field belongs to: "analyses.contacts" and "output.video" for the
// fields of an analysis or an output, else the first key ("inputs", "visual", "notes").
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

}  // namespace

PlaceholderPage::PlaceholderPage(forms::FormContext& context, const QString& id, QWidget* parent)
    : QWidget(parent), m_context(context), m_id(id) {
    setObjectName(QStringLiteral("page:") + id);
    m_layout = new QVBoxLayout(this);
    m_layout->setContentsMargins(6, 6, 6, 6);
    rebuild();
    connect(&m_context, &forms::FormContext::documentChanged, this, [this](const core::DocumentChange& change) {
        if (change.kind == core::DocumentChange::Kind::Replaced && pagePaths() != m_paths) {
            rebuild();
        }
    });
}

QString PlaceholderPage::title() const { return toQString(core::problemPageTitle(m_id.toStdString())); }

std::vector<std::string> PlaceholderPage::pagePaths() const {
    const core::Project& project = m_context.document().project();
    std::vector<std::string> out;
    std::vector<std::string> lists;  // lists of records shown by a list editor: their items are its cells
    const std::string page = m_id.toStdString();
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

void PlaceholderPage::rebuild() {
    delete m_body;
    m_paths = pagePaths();
    m_body = new QWidget(this);
    auto* layout = new QVBoxLayout(m_body);
    layout->setContentsMargins(0, 0, 0, 0);
    auto* note = new QLabel(tr("All settings of this page, one field each (the page's own layout comes with a later "
                               "version)."),
                            m_body);
    note->setWordWrap(true);
    note->setEnabled(false);  // a hint
    layout->addWidget(note);

    const core::Project& project = m_context.document().project();
    std::map<std::string, PageSection*> sections;
    for (const std::string& path : m_paths) {
        const std::string part = sectionPath(path);
        PageSection*& section = sections[part];
        if (!section) {
            const std::size_t dot = part.rfind('.');
            QString id = toQString(dot == std::string::npos ? part : part.substr(dot + 1));
            QString title = lastPart(core::describeField(part, &project));
            if (part == "analyses" || part == "output") {  // the fields directly under analyses / output
                id = QStringLiteral("general");
                title = part == "analyses" ? tr("General") : tr("Locations");
            }
            section = new PageSection(m_id, id, title, m_body);
            layout->addWidget(section);
        }
        forms::FieldWidget* field = forms::createFieldWidget(m_context, path, section->content());
        if (!field) continue;
        QString label = field->label();
        if (label.isEmpty()) label = lastPart(core::describeField(path, &project));
        section->form()->addRow(label, field);
    }
    layout->addStretch(1);
    m_layout->addWidget(m_body);
}

bool PlaceholderPage::revealField(const QString& path) {
    for (auto* section : m_body->findChildren<PageSection*>()) {  // the field may be in a folded section
        const std::string wanted = path.toStdString();
        for (forms::FieldWidget* field : section->findChildren<forms::FieldWidget*>()) {
            if (field->path() == wanted || inside(wanted, field->path())) section->setExpanded(true);
        }
    }
    return problems::revealFieldIn(this, path);
}

}  // namespace mdbridge::app::pages
