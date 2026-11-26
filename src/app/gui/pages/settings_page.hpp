#pragma once
// What the settings pages share (milestone 5 W1-B): the scrolled page, its collapsible sections and
// MDATK's two-column form (gui/widgets.py Form, frame, hint), fields made from defaultFieldSpec with
// the page's labels and tooltips, and revealField().
//
// A page is MDATK's side-bar page: a vertical list of framed sections in a scroll area that never
// scrolls sideways (the side bar has a fixed width). A section is a group box named
// "section:<page>.<section>" with a bold title button that folds it; its form has a label column
// and a widget column, a hint below a field is a dim wrapped label.

#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

#include <QGroupBox>
#include <QString>
#include <QWidget>

#include "gui/forms/field_spec.hpp"
#include "gui/forms/form_context.hpp"

class QGridLayout;
class QLabel;
class QPushButton;
class QScrollArea;
class QToolButton;
class QVBoxLayout;

namespace mdbridge::app::forms {
class FieldWidget;
}

namespace mdbridge::app::pages {

class Section : public QGroupBox {
    Q_OBJECT

public:
    // objectName "section:<page>.<id>"
    Section(const QString& page, const QString& id, const QString& title, QWidget* parent = nullptr);

    QString sectionId() const { return m_id; }
    QString title() const { return m_title; }
    QToolButton* header() const { return m_header; }
    bool isExpanded() const;
    void setExpanded(bool expanded);

    // MDATK's Form: a row "label | widget" (the tooltip on both; the widget keeps its own when it has
    // one), a widget across both columns, a row of widgets side by side. Return the widget / row.
    QWidget* add(const QString& label, QWidget* widget, const QString& tooltip = {});
    QWidget* addWide(QWidget* widget);
    QWidget* addRow(const QString& label, std::initializer_list<QWidget*> widgets, const QString& tooltip = {});
    QWidget* body() const { return m_body; }

private:
    QString m_id;
    QString m_title;
    QToolButton* m_header = nullptr;
    QWidget* m_body = nullptr;
    QGridLayout* m_grid = nullptr;
    int m_row = 0;
};

// MDATK's hint(): dim, wrapped text.
QLabel* hintLabel(const QString& text, QWidget* parent = nullptr);
// A plain label between the widgets of a row ("vs", "bonds").
QLabel* rowLabel(const QString& text, QWidget* parent = nullptr);
// A push button named `name` with `text`, tooltip and accessible name.
QPushButton* pageButton(const QString& name, const QString& text, const QString& tooltip = {},
                        QWidget* parent = nullptr);

class SettingsPage : public QWidget {
    Q_OBJECT

public:
    SettingsPage(forms::FormContext& context, QString id, QString title, QWidget* parent = nullptr);

    // One of core::kProblemPages' ids ("input") and its title ("Input & cleaning").
    Q_INVOKABLE QString pageId() const { return m_id; }
    Q_INVOKABLE QString title() const { return m_title; }

    // Scrolls to the widget of `path` and gives it the focus: a field ("inputs.stop"), a list
    // editor ("groups") or one of its cells / rows ("groups.1.residues", "groups.1"), the section's
    // fold opened. False when the page shows no such field.
    Q_INVOKABLE bool revealField(const QString& path);

    forms::FormContext& context() const { return m_context; }
    const std::vector<Section*>& sections() const { return m_sections; }
    Section* section(const QString& id) const;
    QScrollArea* scrollArea() const { return m_scroll; }

signals:
    // A line for the window's log ("WARNING: ..." for problems), as MDATK's append_log.
    void logMessage(const QString& text);

protected:
    Section* addSection(const QString& id, const QString& title);
    // The field widget of `path`: defaultFieldSpec(model path), changed by `adjust`, with `hint` as
    // tooltip when the spec has none and `label` as its label when given (check box text,
    // accessible name).
    forms::FieldWidget* field(const std::string& path, const QString& label = {}, const QString& hint = {},
                              const std::function<void(forms::FieldSpec&)>& adjust = {});
    // Ends the content (a stretch below the last section).
    void finish();

private:
    forms::FormContext& m_context;
    QString m_id;
    QString m_title;
    QScrollArea* m_scroll = nullptr;
    QWidget* m_content = nullptr;
    QVBoxLayout* m_layout = nullptr;
    std::vector<Section*> m_sections;
};

}  // namespace mdbridge::app::pages
