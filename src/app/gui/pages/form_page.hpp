#pragma once
// The parts a settings page is built from (milestone 5 wave 1): MDATK's two-column form
// (gui/widgets.py Form: narrow wrapped labels, widgets beside them, rows of several widgets, wide
// rows, dim hints), its framed sections as collapsible group boxes, and the scrolled page that
// holds them with field lookup for revealField (the problems list's "go to the field").
//
// Pages derive from FormPage (a QWidget): pageId() is one of core::kProblemPages' ids, title() its
// title; sections are named "section:<page>.<section>".

#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

#include <QGroupBox>
#include <QString>
#include <QWidget>

#include "gui/forms/field_spec.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"

class QGridLayout;
class QLabel;
class QPushButton;
class QScrollArea;
class QToolButton;
class QVBoxLayout;

namespace mdbridge::app::pages {

// MDATK's Form: a label column (wrapped, about 18 characters wide) and a widget column.
class FormGrid : public QWidget {
    Q_OBJECT

public:
    explicit FormGrid(QWidget* parent = nullptr);

    // A row "label | widget" (Form.add); `tooltip` goes on the label and on the widget when it has
    // none (a field widget's hint belongs in its spec, where its flag tooltip keeps it).
    QWidget* add(const QString& label, QWidget* widget, const QString& tooltip = {});
    // A row "label | w1 w2 ..." (Form.add_row): the widgets side by side, not stretched, moving onto
    // further lines when the side bar is too narrow for them (MDATK 0.1.1 flow rows).
    QWidget* addRow(const QString& label, std::initializer_list<QWidget*> widgets, const QString& tooltip = {});
    // A widget across both columns (Form.add_wide).
    QWidget* addWide(QWidget* widget);
    // A dim, wrapped text across both columns (W.hint).
    QLabel* addHint(const QString& text);
    // A row of buttons across both columns (left aligned, wrapping like addRow).
    QWidget* addButtons(std::initializer_list<QWidget*> buttons);

private:
    QLabel* rowLabel(const QString& text, QWidget* buddy, const QString& tooltip);
    QGridLayout* m_grid = nullptr;
    int m_row = 0;
};

// A framed section (MDATK's W.frame) that folds: a header button with the bold title and an
// arrow, the form below it. Named "section:<page>.<id>".
class Section : public QGroupBox {
    Q_OBJECT

public:
    Section(const QString& page, const QString& id, const QString& heading, QWidget* parent = nullptr);

    FormGrid* form() const { return m_form; }
    QToolButton* header() const { return m_header; }
    QString title() const { return m_title; }
    bool isExpanded() const;
    void setExpanded(bool expanded);

private:
    QString m_title;
    QToolButton* m_header = nullptr;
    FormGrid* m_form = nullptr;
};

// A push button with MDATK's tooltip, named `name` (W.button).
QPushButton* makeButton(const QString& text, const QString& name, const QString& tooltip = {},
                        QWidget* parent = nullptr);

// A dim, wrapped label (W.hint).
QLabel* makeHint(const QString& text, QWidget* parent = nullptr);

// A settings page: a scroll area (horizontal scrolling off, as MDATK's W.scrolled) holding the
// sections.
class FormPage : public QWidget {
    Q_OBJECT

public:
    FormPage(forms::FormContext& context, QString id, QString pageTitle, QWidget* parent = nullptr);

    Q_INVOKABLE QString pageId() const { return m_id; }
    Q_INVOKABLE QString title() const { return m_title; }
    // Scrolls to the widget of `path` (a field path; a list editor's cell "visual.representations.1.size";
    // a path no widget shows goes to the nearest field above or beside it: "density.transform" to
    // the density map's fields), opens its section and gives it the focus. False when the path
    // belongs to another page.
    Q_INVOKABLE bool revealField(const QString& path);

    forms::FormContext& context() const { return m_context; }
    QScrollArea* scrollArea() const { return m_scroll; }
    const std::vector<Section*>& sections() const { return m_sections; }
    Section* section(const QString& id) const;

protected:
    // Adds a section at the end of the page.
    Section* addSection(const QString& id, const QString& title);
    // The field widget of `path` with defaultFieldSpec() changed by `adjust` (labels, hints,
    // choices). Throws std::logic_error when the path names no field (a page bug).
    forms::FieldWidget* field(const std::string& path, const std::function<void(forms::FieldSpec&)>& adjust = {});
    // A check box field with MDATK's text (and tooltip).
    forms::FieldWidget* check(const std::string& path, const QString& text, const QString& tooltip = {});

private:
    forms::FormContext& m_context;
    QString m_id;
    QString m_title;
    QScrollArea* m_scroll = nullptr;
    QWidget* m_body = nullptr;
    QVBoxLayout* m_column = nullptr;
    std::vector<Section*> m_sections;
};

}  // namespace mdbridge::app::pages
