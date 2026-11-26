#pragma once
// The field widget for a field path (milestone 5 W0-C): the editor the spec names, or the one the
// field's kind calls for (core::fieldInfo): Integer / Number -> NumberField; Text, optional values,
// comma-separated lists, free-form tables (JSON) and free values -> LineField; Boolean ->
// CheckField; lists of records and number rows -> ListEditor.

#include <string>
#include <vector>

#include "gui/forms/field_spec.hpp"
#include "gui/forms/field_widgets.hpp"

namespace mdbridge::app::forms {

// The editor createFieldWidget() builds for `path` with `spec` (Editor::Auto resolved); Auto when
// the path names nothing or a record (no widget edits a record as a whole).
Editor editorFor(const core::Project& project, const std::string& path, const FieldSpec& spec);

// The widget for `path` with defaultFieldSpec() of its model path, or with `spec`. Null when the
// path names nothing or a record.
FieldWidget* createFieldWidget(FormContext& context, const std::string& path, QWidget* parent = nullptr);
FieldWidget* createFieldWidget(FormContext& context, const std::string& path, const FieldSpec& spec,
                               QWidget* parent = nullptr);

// The fields a complete form of `project` edits: every path of core::fieldPaths() and, before the
// items of a list of records, the list itself (its list editor).
std::vector<std::string> formFieldPaths(const core::Project& project);

}  // namespace mdbridge::app::forms
