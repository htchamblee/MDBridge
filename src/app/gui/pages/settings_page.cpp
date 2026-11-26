#include "gui/pages/settings_page.hpp"

#include <utility>

#include <QAbstractSpinBox>
#include <QFont>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPalette>
#include <QPushButton>
#include <QScrollArea>
#include <QToolButton>
#include <QVBoxLayout>

#include "gui/forms/field_factory.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/problems/problems_panel.hpp"

namespace mdbridge::app::pages {

// ---- Section ---------------------------------------------------------------------------------

Section::Section(const QString& page, const QString& id, const QString& title, QWidget* parent)
    : QGroupBox(parent), m_id(id), m_title(title) {
    setObjectName(QStringLiteral("section:%1.%2").arg(page, id));
    setAccessibleName(title);
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(6, 4, 6, 6);
    outer->setSpacing(2);
    m_header = new QToolButton(this);
    m_header->setObjectName(objectName() + QStringLiteral(":fold"));
    m_header->setText(title);
    m_header->setCheckable(true);
    m_header->setChecked(true);
    m_header->setAutoRaise(true);
    m_header->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    m_header->setArrowType(Qt::DownArrow);
    m_header->setAccessibleName(title);
    m_header->setToolTip(tr("Fold or unfold this section"));
    QFont bold = m_header->font();
    bold.setBold(true);
    m_header->setFont(bold);
    outer->addWidget(m_header);
    m_body = new QWidget(this);
    m_grid = new QGridLayout(m_body);
    m_grid->setContentsMargins(8, 4, 8, 8);  // MDATK's Form: margin 8, column spacing 12, row spacing 6
    m_grid->setHorizontalSpacing(12);
    m_grid->setVerticalSpacing(6);
    m_grid->setColumnStretch(1, 1);
    outer->addWidget(m_body);
    connect(m_header, &QToolButton::toggled, this, &Section::setExpanded);
}

bool Section::isExpanded() const { return m_header->isChecked(); }

void Section::setExpanded(bool expanded) {
    if (m_header->isChecked() != expanded) {
        m_header->setChecked(expanded);  // comes back here through toggled
        return;
    }
    m_header->setArrowType(expanded ? Qt::DownArrow : Qt::RightArrow);
    m_body->setVisible(expanded);
}

QWidget* Section::add(const QString& label, QWidget* widget, const QString& tooltip) {
    auto* text = new QLabel(label, m_body);
    text->setWordWrap(true);
    text->setMaximumWidth(150);  // narrow label column: the side bar has a fixed width (MDATK: 18 chars)
    if (!tooltip.isEmpty()) {
        text->setToolTip(tooltip);
        // field widgets show their own hint (and problems); other widgets take the row's
        if (!qobject_cast<forms::FieldWidget*>(widget) && widget->toolTip().isEmpty()) {
            widget->setToolTip(tooltip);
        }
    }
    if (auto* field = qobject_cast<forms::FieldWidget*>(widget)) {
        text->setBuddy(field->editor());
    }
    widget->setParent(m_body);
    m_grid->addWidget(text, m_row, 0, Qt::AlignLeft | Qt::AlignVCenter);
    m_grid->addWidget(widget, m_row, 1);
    ++m_row;
    return widget;
}

QWidget* Section::addWide(QWidget* widget) {
    widget->setParent(m_body);
    m_grid->addWidget(widget, m_row, 0, 1, 2);
    ++m_row;
    return widget;
}

QWidget* Section::addRow(const QString& label, std::initializer_list<QWidget*> widgets, const QString& tooltip) {
    auto* row = new QWidget(m_body);
    auto* box = new QHBoxLayout(row);
    box->setContentsMargins(0, 0, 0, 0);
    box->setSpacing(6);
    for (QWidget* w : widgets) {
        if (auto* field = qobject_cast<forms::FieldWidget*>(w)) {
            if (auto* spin = qobject_cast<QAbstractSpinBox*>(field->editor())) {
                spin->setMinimumWidth(72);  // MDATK: 5 characters and the arrows
            }
        }
        box->addWidget(w);
    }
    box->addStretch(1);
    add(label, row, tooltip);
    return row;
}

QLabel* hintLabel(const QString& text, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    QPalette palette = label->palette();
    palette.setColor(QPalette::WindowText, palette.color(QPalette::PlaceholderText));  // MDATK's dim-label
    label->setPalette(palette);
    return label;
}

QLabel* rowLabel(const QString& text, QWidget* parent) { return new QLabel(text, parent); }

QPushButton* pageButton(const QString& name, const QString& text, const QString& tooltip, QWidget* parent) {
    auto* button = new QPushButton(text, parent);
    button->setObjectName(name);
    button->setAccessibleName(text);
    if (!tooltip.isEmpty()) {
        button->setToolTip(tooltip);
        button->setAccessibleDescription(tooltip);
    }
    return button;
}

// ---- SettingsPage ----------------------------------------------------------------------------

SettingsPage::SettingsPage(forms::FormContext& context, QString id, QString title, QWidget* parent)
    : QWidget(parent), m_context(context), m_id(std::move(id)), m_title(std::move(title)) {
    setObjectName(QStringLiteral("page:") + m_id);
    setAccessibleName(m_title);
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(0, 0, 0, 0);
    m_scroll = new QScrollArea(this);
    m_scroll->setObjectName(QStringLiteral("page:%1:scroll").arg(m_id));
    m_scroll->setWidgetResizable(true);
    m_scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);  // rows fit the side bar, never scroll sideways
    m_scroll->setFrameShape(QFrame::NoFrame);
    m_content = new QWidget;
    m_layout = new QVBoxLayout(m_content);
    m_layout->setContentsMargins(10, 10, 10, 10);
    m_layout->setSpacing(8);
    m_scroll->setWidget(m_content);
    outer->addWidget(m_scroll);
}

Section* SettingsPage::addSection(const QString& id, const QString& title) {
    auto* s = new Section(m_id, id, title, m_content);
    m_layout->addWidget(s);
    m_sections.push_back(s);
    return s;
}

Section* SettingsPage::section(const QString& id) const {
    for (Section* s : m_sections) {
        if (s->sectionId() == id) {
            return s;
        }
    }
    return nullptr;
}

void SettingsPage::finish() { m_layout->addStretch(1); }

forms::FieldWidget* SettingsPage::field(const std::string& path, const QString& label, const QString& hint,
                                        const std::function<void(forms::FieldSpec&)>& adjust) {
    const auto info = m_context.document().info(path);
    forms::FieldSpec spec = forms::defaultFieldSpec(info ? info->modelPath : path);
    if (!label.isEmpty()) {
        spec.label = label;
    }
    if (spec.hint.isEmpty()) {
        spec.hint = hint;
    }
    if (adjust) {
        adjust(spec);
    }
    forms::FieldWidget* w = forms::createFieldWidget(m_context, path, spec, m_content);
    if (!w) {
        // a page names only fields the model has; a missing one is a programming error
        qFatal("settings page %s: no field widget for '%s'", qPrintable(m_id), path.c_str());
    }
    return w;
}

namespace {

Section* sectionOf(QWidget* w) {
    for (QWidget* p = w; p; p = p->parentWidget()) {
        if (auto* s = qobject_cast<Section*>(p)) {
            return s;
        }
    }
    return nullptr;
}

}  // namespace

bool SettingsPage::revealField(const QString& path) {
    // open the fold of the section that holds the field (or the nearest parent path's widget)
    const std::string key = forms::toStdString(path);
    for (std::string p = key; !p.empty();) {
        QWidget* target = nullptr;
        for (forms::FieldWidget* f : m_content->findChildren<forms::FieldWidget*>()) {
            if (f->path() == p) {
                target = f;
                break;
            }
        }
        if (!target) {
            target = forms::findCell(m_content, p).editor;
        }
        if (target) {
            if (Section* s = sectionOf(target)) {
                s->setExpanded(true);
            }
            break;
        }
        const std::size_t dot = p.rfind('.');
        p = dot == std::string::npos ? std::string() : p.substr(0, dot);
    }
    return problems::revealFieldIn(m_content, path, true);
}

}  // namespace mdbridge::app::pages
