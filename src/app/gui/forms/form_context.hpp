#pragma once
// The bound field widgets of one document (milestone 5 W0-C): the document's changes as a Qt
// signal, its problems (computed once per change, not per widget) and the form errors of text a
// field could not take (AP-10), and the fields that hold typed text not yet stored.
//
// The window keeps one FormContext per open document; the document must outlive it, and the
// context must outlive its widgets.

#include <map>
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

    core::ProjectDocument& m_document;
    core::ProjectDocument::ListenerId m_listener = 0;
    std::vector<core::Problem> m_problems;
    bool m_problemsDirty = true;
    std::map<std::string, core::Problem> m_formErrors;
    std::vector<QPointer<FieldWidget>> m_fields;
    QTimer m_problemTimer;
};

}  // namespace mdbridge::app::forms
