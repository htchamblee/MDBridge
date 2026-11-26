#include "gui/forms/form_context.hpp"

#include <algorithm>
#include <exception>
#include <set>
#include <tuple>
#include <utility>

#include "gui/forms/field_widgets.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_file.hpp"

namespace mdbridge::app::forms {

QString toQString(std::string_view text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }

std::string toStdString(const QString& text) {
    const QByteArray bytes = text.toUtf8();
    return std::string(bytes.constData(), static_cast<std::size_t>(bytes.size()));
}

FormContext::FormContext(core::ProjectDocument& document, QObject* parent) : QObject(parent), m_document(document) {
    m_problemTimer.setSingleShot(true);
    m_problemTimer.setInterval(0);
    connect(&m_problemTimer, &QTimer::timeout, this, &FormContext::refreshProblems);
    m_listener = m_document.addListener([this](const core::DocumentChange& change) { onChange(change); });
}

FormContext::~FormContext() { m_document.removeListener(m_listener); }

void FormContext::onChange(const core::DocumentChange& change) {
    m_problemsDirty = true;
    // a field that shows the document's value again has no form error any more (undo, open)
    for (auto it = m_formErrors.begin(); it != m_formErrors.end();) {
        it = change.affects(it->first) ? m_formErrors.erase(it) : std::next(it);
    }
    emit documentChanged(change);
    // one check per change; a drag (an undo group) is checked once it pauses
    m_problemTimer.start(m_document.inGroup() ? kCheckDelayMs : 0);
}

void FormContext::scheduleCheck() { m_problemTimer.start(kCheckDelayMs); }

std::vector<core::Problem> FormContext::computeProblems() const {
    // the typed text not yet stored (MDATK 0.1.1 checks the project as the side bar shows it)
    std::vector<std::pair<std::string, std::string>> typed;
    for (const QPointer<FieldWidget>& field : m_fields) {
        if (field) {
            if (const std::optional<QString> text = field->pendingText()) {
                typed.emplace_back(field->path(), toStdString(*text));
            }
        }
    }
    if (typed.empty()) {
        return m_document.problems();
    }
    core::ProjectDocument scratch;  // no file: nothing is read from disk
    scratch.replace(m_document.project());
    std::vector<core::Problem> errors;
    for (const auto& [path, text] : typed) {
        try {
            scratch.setText(path, text);
        } catch (const core::FieldError& error) {
            errors.push_back(error.problem(&scratch.project()));
        } catch (const std::exception& error) {
            errors.push_back(core::locatedProblem(path, error.what(), &scratch.project()));
        }
    }
    core::Project resolved = scratch.project();
    try {
        core::resolvePaths(resolved, m_document.folder());
    } catch (const core::ProjectError&) {
        resolved = scratch.project();  // as ProjectDocument::problems: the project as written
    }
    std::vector<core::Problem> found = core::projectProblems(resolved);
    errors.insert(errors.end(), std::make_move_iterator(found.begin()), std::make_move_iterator(found.end()));
    return errors;
}

const std::vector<core::Problem>& FormContext::problems() {
    // inside a drag the last check stands until the drag pauses (no check per mouse move)
    if (m_problemsDirty && (!m_problemsKnown || !m_document.inGroup())) {
        m_problemsDirty = false;
        m_problemsKnown = true;
        try {
            m_problems = computeProblems();
        } catch (const std::exception&) {
            m_problems.clear();  // a project the checks cannot read has no located problems
        }
    }
    return m_problems;
}

void FormContext::refreshProblems() {
    m_problemTimer.stop();
    m_problemsDirty = true;
    m_problemsKnown = false;
    problems();
    emit problemsChanged();
}

std::vector<core::Problem> FormContext::problemsAt(std::string_view path, bool inside) {
    const auto matches = [&](const std::string& field) {
        if (field == path) {
            return true;
        }
        return inside && field.size() > path.size() && field.compare(0, path.size(), path) == 0 &&
               (path.empty() || field[path.size()] == '.');
    };
    std::vector<core::Problem> out;
    std::set<std::tuple<std::string, std::string, std::string>> seen;
    const auto add = [&](std::vector<core::Problem>& to, const core::Problem& p) {
        if (seen.emplace(p.field, p.message, p.level).second) {
            to.push_back(p);
        }
    };
    for (const auto& [field, problem] : m_formErrors) {
        if (matches(field)) {
            add(out, problem);
        }
    }
    std::vector<core::Problem> warnings;
    for (const core::Problem& p : problems()) {
        if (matches(p.field)) {
            add(p.level == "warning" ? warnings : out, p);
        }
    }
    out.insert(out.end(), warnings.begin(), warnings.end());
    return out;
}

std::vector<core::Problem> FormContext::allProblems() {
    std::vector<core::Problem> out;
    std::set<std::pair<std::string, std::string>> seen;
    std::set<std::string> fields;
    const auto add = [&](const core::Problem& p) {
        if (seen.emplace(p.field, p.message).second) {
            out.push_back(p);
            fields.insert(p.field);
        }
    };
    for (const auto& [field, problem] : m_formErrors) {
        add(problem);
    }
    for (const core::Problem& p : problems()) {
        add(p);
    }
    // the widgets' own problems show where a field has no other (as the field's flag does)
    std::vector<core::Problem> local;
    for (const auto& [field, problem] : m_localProblems) {
        if (!fields.count(field)) {
            local.push_back(problem);
        }
    }
    for (const core::Problem& p : local) {
        add(p);
    }
    return out;
}

void FormContext::setLocalProblem(const std::string& path, const core::Problem* problem) {
    const auto it = m_localProblems.find(path);
    if (!problem) {
        if (it == m_localProblems.end()) {
            return;
        }
        m_localProblems.erase(it);
    } else if (it != m_localProblems.end() && it->second == *problem) {
        return;
    } else {
        m_localProblems[path] = *problem;
    }
    emit problemsChanged();
}

void FormContext::setFormError(const std::string& path, const core::Problem* problem) {
    const auto it = m_formErrors.find(path);
    if (!problem) {
        if (it == m_formErrors.end()) {
            return;
        }
        m_formErrors.erase(it);
    } else if (it != m_formErrors.end() && it->second == *problem) {
        return;
    } else {
        m_formErrors[path] = *problem;
    }
    emit problemsChanged();
}

void FormContext::commitPendingEdits() {
    const auto fields = m_fields;  // a commit may rebuild widgets
    for (const QPointer<FieldWidget>& field : fields) {
        if (field) {
            field->commit();
        }
    }
}

void FormContext::addField(FieldWidget* field) { m_fields.emplace_back(field); }

void FormContext::removeField(FieldWidget* field) {
    // a field that is gone shows no problem of its own (told once the event loop runs: the field
    // is being destroyed, and a form being closed removes many)
    if (field && m_localProblems.erase(field->path()) && !m_problemTimer.isActive()) {
        m_problemTimer.start(0);
    }
    std::erase_if(m_fields, [field](const QPointer<FieldWidget>& f) { return f.isNull() || f.data() == field; });
}

}  // namespace mdbridge::app::forms
