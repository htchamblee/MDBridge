#pragma once
// Field widgets bound to a ProjectDocument by field path (milestone 5 W0-C; MDATK's gui/widgets.py
// Binder and factories, with the fixes of docs/specs/m5-plan.md AP-01 .. AP-05 and AP-15).
//
// Each FieldWidget edits one field path. It shows the document's value (and again after every
// change of it: undo, redo, open), never writes while showing, and stores what the user enters as
// one undo step: text when editing finishes (Return, leaving the field, commitPendingEdits()), a
// pick at once. Text the field cannot take stays shown and is flagged with its located message
// (AP-10); the document keeps its value. Problems of the document for the field flag it too
// (tinted, an icon in line edits, the message as tooltip above the field's hint, MDATK 0.1.1's
// set_flag).
//
// The editor inside (QLineEdit, the spin box, QComboBox, QCheckBox, QPlainTextEdit, QTableView) has
// the field path as objectName, so tests and the GUI harness find it; the FieldWidget itself is
// "<path>:field", its buttons "<path>:browse", "<path>:pick", "<path>:add" ... Value widgets ignore
// the mouse wheel unless they have the keyboard focus (MDATK 0.0.4: scrolling a page never changes a
// value), and every control is reachable with the keyboard and has an accessible name.

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <QAbstractSpinBox>
#include <QCheckBox>
#include <QColor>
#include <QIcon>
#include <QComboBox>
#include <QPlainTextEdit>
#include <QPointer>
#include <QString>
#include <QStringList>
#include <QWidget>

#include "gui/forms/field_spec.hpp"
#include "gui/forms/form_context.hpp"
#include "mdbridge/core/project_document.hpp"

class QHBoxLayout;
class QLineEdit;
class QToolButton;

namespace mdbridge::app::forms {

// ---- helpers ---------------------------------------------------------------------------------

// The widget ignores wheel events unless it has the keyboard focus (the event goes on to the
// scrolled page); it no longer takes the focus from the wheel. Used on spin boxes, combos and
// sliders (MDATK's no_wheel); hasWheelGuard() tells whether a widget has it.
void guardWheel(QWidget* widget);
bool hasWheelGuard(const QWidget* widget);

// Shows `message` as a problem of `widget` ("error" or "warning"; "" clears it): the widget's base
// colour tinted, an icon at the end of a line edit, and the message as tooltip and accessible
// description above `hint`, which is the tooltip again once the flag is cleared.
void setFlag(QWidget* widget, const QString& message, const QString& level, const QString& hint);
// The (message, level) a widget shows ({"", ""} when none).
std::pair<QString, QString> flagOf(const QWidget* widget);

// The colour Matplotlib reads from `text` (invalid when none) and a swatch of it (crossed out
// when invalid; transparency over a checkerboard).
QColor matplotlibColor(const QString& text);
QIcon colorSwatch(const QColor& colour);

// The dialogs the widgets open. Tests replace them (no native dialog can be driven offscreen).
struct FileDialogRequest {
    PathMode mode = PathMode::OpenFile;
    QString title;
    QString start;        // the folder or file the dialog starts at
    QStringList filters;  // Qt name filters (nameFilters())
};
using FileDialogFunction = std::function<QStringList(QWidget* parent, const FileDialogRequest& request)>;
using ColorDialogFunction = std::function<QColor(QWidget* parent, const QColor& initial, const QString& title)>;
// Native dialogs by default (owner decision D3); returns the previous function.
FileDialogFunction setFileDialogFunction(FileDialogFunction function);
ColorDialogFunction setColorDialogFunction(ColorDialogFunction function);
QStringList runFileDialog(QWidget* parent, const FileDialogRequest& request);
QColor runColorDialog(QWidget* parent, const QColor& initial, const QString& title);

// The text a path picked in a file dialog is stored as (owner decision, in the spirit of D8): when
// the document has a file and `picked` lies inside its folder (the folder itself included), the
// path relative to that folder with '/' separators, as a user would write it ("inputs/a.parm7",
// "."; "./~x" for a name starting with '~', which would read as the home folder); else `picked`
// with the platform's separators (outside the folder, another drive or share, a new project).
QString pickedPathText(const core::ProjectDocument& document, const QString& picked);

// ---- the base --------------------------------------------------------------------------------

class FieldWidget : public QWidget {
    Q_OBJECT

public:
    FieldWidget(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    ~FieldWidget() override;

    const std::string& path() const { return m_path; }
    const FieldSpec& spec() const { return m_spec; }
    FormContext& context() const { return m_context; }
    core::ProjectDocument& document() const { return m_context.document(); }
    // The widget named after the path.
    QWidget* editor() const { return m_editor; }

    // Shows the document's value (no edit). Called after each change that affects the path.
    virtual void refresh() = 0;
    // Stores text typed but not yet stored (nothing when there is none).
    virtual void commit() {}
    // The text the field shows for its value (what a commit of it would store).
    virtual QString displayedText() const = 0;

    // The flag shown: (message, level).
    std::pair<QString, QString> flag() const;
    // The accessible name / label (spec.label, else the last part of the field's location).
    QString label() const;

signals:
    // The user's input was stored (one undo step) or refused (a form error).
    void committed(bool stored);

protected:
    // Names `editor` after the path, gives it the label, hint and placeholder; the flag goes on
    // `flagTarget` (default: the editor).
    void setEditor(QWidget* editor, QWidget* flagTarget = nullptr);
    // Stores `text` (setText, one step; nothing when it is the text shown). False when refused:
    // the form error is flagged.
    bool commitText(const QString& text, const std::string& path = {});
    // The document's text of the field ("" when the path names nothing now).
    QString documentText() const;
    bool fieldExists() const;
    // A problem the widget itself sees (a value outside the range hint, an unknown choice, a
    // colour Matplotlib cannot read); shown when the document has none for the field.
    void setLocalProblem(const QString& message, const QString& level = QStringLiteral("warning"));
    void updateFlag();
    // True while refresh() writes into the editor (its signals are no user input then).
    bool refreshing() const { return m_refreshing > 0; }

    struct RefreshScope {
        explicit RefreshScope(FieldWidget& w) : widget(w) { ++widget.m_refreshing; }
        ~RefreshScope() { --widget.m_refreshing; }
        RefreshScope(const RefreshScope&) = delete;
        RefreshScope& operator=(const RefreshScope&) = delete;
        FieldWidget& widget;
    };

    virtual void onDocumentChanged(const core::DocumentChange& change);

private:
    FormContext& m_context;
    std::string m_path;
    FieldSpec m_spec;
    QWidget* m_editor = nullptr;
    QWidget* m_flagTarget = nullptr;
    QString m_localMessage;
    QString m_localLevel;
    int m_refreshing = 0;
};

// ---- number ----------------------------------------------------------------------------------

// A spin box whose text is any text: no range, no clamping (AP-01), the arrows step by spec.step
// within the range hint.
class NumberBox : public QAbstractSpinBox {
    Q_OBJECT

public:
    explicit NumberBox(QWidget* parent = nullptr);
    QLineEdit* edit() const { return lineEdit(); }
    void stepBy(int steps) override;
    std::function<void(int)> onStep;
    std::function<StepEnabled()> onStepEnabled;
    StepEnabled stepEnabled() const override;
};

// Integer and Number fields. The text is the project writer's (shortest round trip, AP-02); a value
// outside the range hint is kept and flagged as a warning.
class NumberField : public FieldWidget {
    Q_OBJECT

public:
    NumberField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    NumberBox* box() const { return m_box; }
    void refresh() override;
    void commit() override;
    QString displayedText() const override;

private:
    void step(int steps);
    bool isInteger() const;
    std::optional<double> currentNumber() const;
    NumberBox* m_box = nullptr;
};

// ---- text ------------------------------------------------------------------------------------

// One-line text: Text and OptionalText, OptionalInteger / OptionalNumber ("" unset, shown with the
// placeholder), comma-separated lists (TextList "a, b", NumberList), JSON of free-form tables
// (panel options) and free values. Stored when editing finishes; Esc shows the stored value again.
class LineField : public FieldWidget {
    Q_OBJECT

public:
    LineField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QLineEdit* lineEdit() const { return m_edit; }
    void refresh() override;
    void commit() override;
    QString displayedText() const override;

protected:
    virtual void textShown(const QString& text);  // after refresh and while typing

    QLineEdit* m_edit = nullptr;
    QHBoxLayout* m_row = nullptr;
};

// A path: the text as written (relative, "~", absolute, AP-06) and a Browse button with MDATK's
// filters; the file picked is stored at once (native separators).
class PathField : public LineField {
    Q_OBJECT

public:
    PathField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QToolButton* browseButton() const { return m_browse; }
    // Opens the dialog (the button's action).
    void browse();

private:
    QToolButton* m_browse = nullptr;
};

// A Matplotlib colour as text, kept as typed ("tab:orange", "#ff000080", "gold", "0.5"; AP-03),
// with a swatch of the colour and a picker that writes "#rrggbb". Text Matplotlib cannot read is
// flagged with its message (RN-12).
class ColorField : public LineField {
    Q_OBJECT

public:
    ColorField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QToolButton* pickButton() const { return m_pick; }
    // Opens the colour dialog (the swatch's action).
    void pick();
    // The swatch's colour (invalid when the text is no colour).
    QColor swatchColor() const { return m_swatch; }

protected:
    void textShown(const QString& text) override;

private:
    QToolButton* m_pick = nullptr;
    QColor m_swatch;
};

// Multi-line text (notes, extra cpptraj commands): stored when the field loses the focus or on
// commitPendingEdits().
class MultiLineField : public FieldWidget {
    Q_OBJECT

public:
    MultiLineField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QPlainTextEdit* textEdit() const { return m_edit; }
    void refresh() override;
    void commit() override;
    QString displayedText() const override;

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private:
    QPlainTextEdit* m_edit = nullptr;
    bool m_dirty = false;
};

// ---- choices ---------------------------------------------------------------------------------

// Fixed choices (spec.choices: value, label). A value not among them is kept and shown as
// "<value> (unknown)" and flagged (AP-05); picking an item stores it at once.
class ComboField : public FieldWidget {
    Q_OBJECT

public:
    ComboField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QComboBox* comboBox() const { return m_combo; }
    void setChoices(const std::vector<Choice>& choices);
    void refresh() override;
    QString displayedText() const override;

private:
    QComboBox* m_combo = nullptr;
    std::vector<Choice> m_choices;
};

// Free text with suggestions (group names, chain palettes, codecs): stored when editing finishes
// or an item is picked; the text is never replaced by a suggestion.
class EditableComboField : public FieldWidget {
    Q_OBJECT

public:
    EditableComboField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QComboBox* comboBox() const { return m_combo; }
    // Replaces the suggestions, keeping the text (set_combo_entry_items).
    void setChoices(const QStringList& items);
    void refresh() override;
    void commit() override;
    QString displayedText() const override;

private:
    QComboBox* m_combo = nullptr;
};

// ---- check -----------------------------------------------------------------------------------

// Boolean (and OptionalBoolean: partly checked while unset; the user's clicks set true or false).
class CheckField : public FieldWidget {
    Q_OBJECT

public:
    CheckField(FormContext& context, std::string path, FieldSpec spec, QWidget* parent = nullptr);
    QCheckBox* checkBox() const { return m_check; }
    void refresh() override;
    QString displayedText() const override;

private:
    QCheckBox* m_check = nullptr;
};

}  // namespace mdbridge::app::forms
