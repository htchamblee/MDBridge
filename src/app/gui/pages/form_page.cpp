#include "gui/pages/form_page.hpp"

#include <stdexcept>
#include <string_view>
#include <utility>

#include <QAbstractSpinBox>
#include <QFont>
#include <QGridLayout>
#include <QLabel>
#include <QPalette>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QToolButton>
#include <QVBoxLayout>

#include "gui/forms/field_factory.hpp"
#include "gui/forms/flow_layout.hpp"
#include "gui/problems/problems_panel.hpp"
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

namespace {

constexpr int kLabelWidth = 130;   // MDATK: set_max_width_chars(18) in the 520 px side bar
constexpr int kNumberWidth = 84;   // a spin box in a row of several (MDATK's spins show 5-7 digits)

}  // namespace

// ---- FormGrid --------------------------------------------------------------------------------

FormGrid::FormGrid(QWidget* parent) : QWidget(parent) {
    m_grid = new QGridLayout(this);
    m_grid->setContentsMargins(8, 8, 8, 8);
    m_grid->setHorizontalSpacing(12);
    m_grid->setVerticalSpacing(6);
    m_grid->setColumnStretch(1, 1);
}

QLabel* FormGrid::rowLabel(const QString& text, QWidget* buddy, const QString& tooltip) {
    auto* label = new QLabel(text, this);
    label->setWordWrap(true);
    label->setFixedWidth(kLabelWidth);
    label->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    if (!tooltip.isEmpty()) {
        label->setToolTip(tooltip);
    }
    if (buddy) {
        QWidget* target = buddy;
        if (auto* f = qobject_cast<forms::FieldWidget*>(buddy); f && f->editor()) {
            target = f->editor();
        }
        label->setBuddy(target);
    }
    return label;
}

QWidget* FormGrid::add(const QString& label, QWidget* widget, const QString& tooltip) {
    widget->setParent(this);
    if (!tooltip.isEmpty() && !qobject_cast<forms::FieldWidget*>(widget) && widget->toolTip().isEmpty()) {
        widget->setToolTip(tooltip);
    }
    m_grid->addWidget(rowLabel(label, widget, tooltip), m_row, 0);
    m_grid->addWidget(widget, m_row, 1);
    ++m_row;
    return widget;
}

QWidget* FormGrid::addRow(const QString& label, std::initializer_list<QWidget*> widgets, const QString& tooltip) {
    auto* row = new QWidget(this);
    auto* box = new forms::FlowLayout(row);  // wraps when the side bar is too narrow for the row
    box->setContentsMargins(0, 0, 0, 0);
    for (QWidget* w : widgets) {
        if (auto* f = qobject_cast<forms::FieldWidget*>(w); f && qobject_cast<QAbstractSpinBox*>(f->editor())) {
            f->setMinimumWidth(kNumberWidth);
        }
        box->addWidget(w);
        if (!tooltip.isEmpty() && !qobject_cast<forms::FieldWidget*>(w) && w->toolTip().isEmpty()) {
            w->setToolTip(tooltip);
        }
    }
    m_grid->addWidget(rowLabel(label, widgets.size() ? *widgets.begin() : nullptr, tooltip), m_row, 0);
    m_grid->addWidget(row, m_row, 1);
    ++m_row;
    return row;
}

QWidget* FormGrid::addWide(QWidget* widget) {
    widget->setParent(this);
    m_grid->addWidget(widget, m_row, 0, 1, 2);
    ++m_row;
    return widget;
}

QLabel* FormGrid::addHint(const QString& text) {
    QLabel* hint = makeHint(text, this);
    addWide(hint);
    return hint;
}

QWidget* FormGrid::addButtons(std::initializer_list<QWidget*> buttons) {
    auto* row = new QWidget(this);
    auto* box = new forms::FlowLayout(row);
    box->setContentsMargins(0, 0, 0, 0);
    for (QWidget* b : buttons) {
        box->addWidget(b);
    }
    return addWide(row);
}

QLabel* makeHint(const QString& text, QWidget* parent) {
    auto* hint = new QLabel(text, parent);
    hint->setWordWrap(true);
    hint->setForegroundRole(QPalette::PlaceholderText);  // GTK's dim-label
    return hint;
}

QPushButton* makeButton(const QString& text, const QString& name, const QString& tooltip, QWidget* parent) {
    auto* button = new QPushButton(text, parent);
    button->setObjectName(name);
    if (!tooltip.isEmpty()) {
        button->setToolTip(tooltip);
    }
    return button;
}

// ---- Section ---------------------------------------------------------------------------------

Section::Section(const QString& page, const QString& id, const QString& heading, QWidget* parent)
    : QGroupBox(parent), m_title(heading) {
    setObjectName(QStringLiteral("section:%1.%2").arg(page, id));
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(4, 4, 4, 4);
    layout->setSpacing(0);
    m_header = new QToolButton(this);
    m_header->setObjectName(objectName() + QStringLiteral(":header"));
    m_header->setText(heading);
    m_header->setCheckable(true);
    m_header->setChecked(true);
    m_header->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    m_header->setArrowType(Qt::DownArrow);
    m_header->setAutoRaise(true);
    m_header->setAccessibleName(heading);
    QFont bold = m_header->font();
    bold.setBold(true);
    m_header->setFont(bold);
    m_form = new FormGrid(this);
    layout->addWidget(m_header, 0, Qt::AlignLeft);
    layout->addWidget(m_form);
    connect(m_header, &QToolButton::toggled, this, &Section::setExpanded);
}

bool Section::isExpanded() const { return m_header->isChecked(); }

void Section::setExpanded(bool expanded) {
    if (m_header->isChecked() != expanded) {
        m_header->setChecked(expanded);  // toggled() calls this again
        return;
    }
    m_header->setArrowType(expanded ? Qt::DownArrow : Qt::RightArrow);
    m_form->setVisible(expanded);
}

// ---- FormPage --------------------------------------------------------------------------------

FormPage::FormPage(forms::FormContext& context, QString id, QString pageTitle, QWidget* parent)
    : QWidget(parent), m_context(context), m_id(std::move(id)), m_title(std::move(pageTitle)) {
    setObjectName(QStringLiteral("page:") + m_id);
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_scroll = new QScrollArea(this);
    m_scroll->setObjectName(objectName() + QStringLiteral(":scroll"));
    m_scroll->setWidgetResizable(true);
    m_scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    m_scroll->setFrameShape(QFrame::NoFrame);
    m_body = new QWidget;
    m_body->setObjectName(objectName() + QStringLiteral(":body"));
    m_column = new QVBoxLayout(m_body);
    m_column->setContentsMargins(10, 10, 10, 10);
    m_column->setSpacing(8);
    m_column->addStretch(1);
    m_scroll->setWidget(m_body);
    layout->addWidget(m_scroll);
}

Section* FormPage::addSection(const QString& id, const QString& sectionTitle) {
    auto* s = new Section(m_id, id, sectionTitle, m_body);
    m_column->insertWidget(m_column->count() - 1, s);  // before the stretch
    m_sections.push_back(s);
    return s;
}

Section* FormPage::section(const QString& id) const {
    const QString name = QStringLiteral("section:%1.%2").arg(m_id, id);
    for (Section* s : m_sections) {
        if (s->objectName() == name) {
            return s;
        }
    }
    return nullptr;
}

forms::FieldWidget* FormPage::field(const std::string& path, const std::function<void(forms::FieldSpec&)>& adjust) {
    const auto info = core::fieldInfo(m_context.document().project(), path);
    if (!info) {
        throw std::logic_error("no field '" + path + "' for the " + m_id.toStdString() + " page");
    }
    forms::FieldSpec spec = forms::defaultFieldSpec(info->modelPath);
    if (adjust) {
        adjust(spec);
    }
    forms::FieldWidget* w = forms::createFieldWidget(m_context, path, spec, m_body);
    if (!w) {
        throw std::logic_error("no field widget for '" + path + "'");
    }
    return w;
}

forms::FieldWidget* FormPage::check(const std::string& path, const QString& text, const QString& tooltip) {
    return field(path, [&](forms::FieldSpec& s) {
        s.label = text;
        if (!tooltip.isEmpty()) {
            s.hint = tooltip;
        }
    });
}

bool FormPage::revealField(const QString& path) {
    const std::string_view page = core::problemPage(path.toStdString());
    if (QString::fromUtf8(page.data(), static_cast<qsizetype>(page.size())) != m_id) {
        return false;
    }
    const auto showField = [this](const QString& field) {
        // a folded section that holds the field opens first (hidden widgets take no focus)
        for (Section* s : m_sections) {
            if (!s->isExpanded() && problems::revealFieldIn(s->form(), field, false)) {
                s->setExpanded(true);
            }
        }
        return problems::revealFieldIn(m_body, field, true);
    };
    if (showField(path)) {
        return true;
    }
    // no widget edits the path or a parent of it (density.transform: the placement): the first
    // field of the nearest parent ("density.map_path")
    for (QString level = path; level.contains(QLatin1Char('.'));) {
        level = level.left(level.lastIndexOf(QLatin1Char('.')));
        const QString prefix = level + QLatin1Char('.');
        for (const QWidget* w : m_body->findChildren<QWidget*>()) {
            if (w->objectName().startsWith(prefix) && !w->objectName().contains(QLatin1Char(':')) &&
                showField(w->objectName())) {
                return true;
            }
        }
    }
    return false;
}

}  // namespace mdbridge::app::pages
