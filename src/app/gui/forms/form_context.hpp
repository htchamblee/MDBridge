#pragma once
// The bound field widgets of one document (milestone 5 W0-C): the document's changes as a Qt
// signal, its problems (computed once per change, not per widget) and the form errors of text a
// field could not take (AP-10), and the fields that hold typed text not yet stored.
//
// Checking as you type (W1-D, MDATK 0.1.1's schedule_check / check_problems): the problems are
// those of the project as the fields show it, typed text not yet stored included (text a field
// could not take is a located error there, as MDATK 0.1.1's store_fields errors). Typing schedules
// one check kCheckDelayMs after the last key; a change of the document is checked at once (one
// check per change), except inside an undo group (a drag), which is checked kCheckDelayMs after
// its last change. allProblems() is what the problems list and the page badges show: these, and
// the problems the field widgets see themselves (a value outside the range hint, an unknown
// choice, a colour Matplotlib cannot read) where a field has no other.
//
// The window keeps one FormContext per open document; the document must outlive it, and the
// context must outlive its widgets.

#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <QObject>
#include <QPointer>
#include <QString>
#include <QTimer>

#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::forms {

class FieldWidget;

QString toQString(std::string_view text);
std::string toStdString(const QString& text);

class FormContext : public QObject {
    Q_OBJECT

public:
    explicit FormContext(core::ProjectDocument& document, QObject* parent = nullptr);
    ~FormContext() override;
    FormContext(const FormContext&) = delete;
    FormContext& operator=(const FormContext&) = delete;

    core::ProjectDocument& document() const { return m_document; }

    // The document's problems (ProjectDocument::problems), recomputed at most once per change.
    const std::vector<core::Problem>& problems();
    // The problems of `path` (also those of fields inside it with `inside`): the form error typed
    // there first, then the document's, errors before warnings.
    std::vector<core::Problem> problemsAt(std::string_view path, bool inside = false);
    // Recomputes the problems now (they are recomputed after a change once the event loop runs).
    void refreshProblems();

    // MDATK 0.1.1's CHECK_DELAY_MS: typed text is checked this long after the last key.
    static constexpr int kCheckDelayMs = 400;
    // A field's text was typed: check (again) kCheckDelayMs from now (one check for a burst).
    void scheduleCheck();
    // A check is due (typing or a drag scheduled one that has not run yet).
    bool checkPending() const { return m_problemTimer.isActive(); }

    // A problem a field widget sees itself (the range hint, an unknown choice, a colour), located;
    // cleared with null.
    void setLocalProblem(const std::string& path, const core::Problem* problem);
    const std::map<std::string, core::Problem>& localProblems() const { return m_localProblems; }

    // Everything the problems list shows, in MDATK 0.1.1's order: the form errors (text a field
    // could not take), the problems of the project as the fields show it (errors in validate()'s
    // order, then warnings), then the field widgets' own problems of fields that have no other;
    // each (field, message) once.
    std::vector<core::Problem> allProblems();

    // Text a field could not take ("'x' is not a whole number"), located; cleared with an empty
    // optional (the field took a value, or shows the document's again).
    void setFormError(const std::string& path, const core::Problem* problem);
    const std::map<std::string, core::Problem>& formErrors() const { return m_formErrors; }

    // Stores the text typed but not yet committed in every field (before a save, a run or closing).
    void commitPendingEdits();

    void addField(FieldWidget* field);
    void removeField(FieldWidget* field);

signals:
    // After every change of the document (undo, redo, open and replace included).
    void documentChanged(const mdbridge::core::DocumentChange& change);
    // The problems or form errors changed (or may have).
    void problemsChanged();

private:
    void onChange(const core::DocumentChange& change);
    std::vector<core::Problem> computeProblems() const;

    core::ProjectDocument& m_document;
    core::ProjectDocument::ListenerId m_listener = 0;
    std::vector<core::Problem> m_problems;
    bool m_problemsDirty = true;
    bool m_problemsKnown = false;
    std::map<std::string, core::Problem> m_formErrors;
    std::map<std::string, core::Problem> m_localProblems;
    std::vector<QPointer<FieldWidget>> m_fields;
    QTimer m_problemTimer;
};

}  // namespace mdbridge::app::forms
