#include "gui/forms/form_context.hpp"

#include <algorithm>
#include <exception>

#include "gui/forms/field_widgets.hpp"

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
    m_problemTimer.start();
}

const std::vector<core::Problem>& FormContext::problems() {
    if (m_problemsDirty) {
        m_problemsDirty = false;
        try {
            m_problems = m_document.problems();
        } catch (const std::exception&) {
            m_problems.clear();  // a project the checks cannot read has no located problems
        }
    }
    return m_problems;
}

void FormContext::refreshProblems() {
    m_problemTimer.stop();
    m_problemsDirty = true;
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
    for (const auto& [field, problem] : m_formErrors) {
        if (matches(field)) {
            out.push_back(problem);
        }
    }
    std::vector<core::Problem> warnings;
    for (const core::Problem& p : problems()) {
        if (matches(p.field)) {
            (p.level == "warning" ? warnings : out).push_back(p);
        }
    }
    out.insert(out.end(), warnings.begin(), warnings.end());
    return out;
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
    std::erase_if(m_fields, [field](const QPointer<FieldWidget>& f) { return f.isNull() || f.data() == field; });
}

}  // namespace mdbridge::app::forms
