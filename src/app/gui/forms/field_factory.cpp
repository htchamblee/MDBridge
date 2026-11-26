#include "gui/forms/field_factory.hpp"

#include <algorithm>
#include <optional>

#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::app::forms {

Editor editorFor(const core::Project& project, const std::string& path, const FieldSpec& spec) {
    const std::optional<core::FieldInfo> info = core::fieldInfo(project, path);
    if (!info || info->kind == core::FieldKind::Record) {
        return Editor::Auto;
    }
    if (spec.editor != Editor::Auto) {
        return spec.editor;
    }
    switch (info->kind) {
    case core::FieldKind::Integer:
    case core::FieldKind::Number:
        return Editor::Number;
    case core::FieldKind::Boolean:
    case core::FieldKind::OptionalBoolean:
        return Editor::Check;
    case core::FieldKind::Records:
    case core::FieldKind::NumberRows:
        return Editor::Table;
    case core::FieldKind::Text:
    case core::FieldKind::OptionalText:
    case core::FieldKind::OptionalInteger:
    case core::FieldKind::OptionalNumber:
    case core::FieldKind::TextList:
    case core::FieldKind::NumberList:
    case core::FieldKind::FreeTable:
    case core::FieldKind::FreeValue:
    case core::FieldKind::Record:
        break;
    }
    return Editor::Line;
}

FieldWidget* createFieldWidget(FormContext& context, const std::string& path, QWidget* parent) {
    const std::optional<core::FieldInfo> info = core::fieldInfo(context.document().project(), path);
    if (!info) {
        return nullptr;
    }
    return createFieldWidget(context, path, defaultFieldSpec(info->modelPath), parent);
}

FieldWidget* createFieldWidget(FormContext& context, const std::string& path, const FieldSpec& spec, QWidget* parent) {
    switch (editorFor(context.document().project(), path, spec)) {
    case Editor::Auto:
        return nullptr;
    case Editor::Number:
        return new NumberField(context, path, spec, parent);
    case Editor::Line:
        return new LineField(context, path, spec, parent);
    case Editor::Path:
        return new PathField(context, path, spec, parent);
    case Editor::Combo:
        return new ComboField(context, path, spec, parent);
    case Editor::EditableCombo:
        return new EditableComboField(context, path, spec, parent);
    case Editor::Color:
        return new ColorField(context, path, spec, parent);
    case Editor::Check:
        return new CheckField(context, path, spec, parent);
    case Editor::MultiLine:
        return new MultiLineField(context, path, spec, parent);
    case Editor::Table:
        return new ListEditor(context, path, spec, parent);
    }
    return nullptr;
}

namespace {

// The lists of records in the project's plain form, empty ones too (fieldPaths names only items).
void recordLists(const core::Project& project, const core::Value& value, const std::string& path,
                 std::vector<std::string>& out) {
    if (const core::Table* table = value.get<core::Table>()) {
        for (const auto& [key, entry] : *table) {
            recordLists(project, entry, path.empty() ? key : path + "." + key, out);
        }
    } else if (const core::Array* array = value.get<core::Array>()) {
        const auto info = core::fieldInfo(project, path);
        if (!info || info->kind != core::FieldKind::Records) {
            return;
        }
        out.push_back(path);
        for (std::size_t i = 0; i < array->size(); ++i) {
            recordLists(project, (*array)[i], path + "." + std::to_string(i), out);
        }
    }
}

}  // namespace

std::vector<std::string> formFieldPaths(const core::Project& project) {
    std::vector<std::string> out;
    for (const std::string& path : core::fieldPaths(project)) {
        // the lists of records above this path, each once
        for (std::size_t dot = path.find('.'); dot != std::string::npos; dot = path.find('.', dot + 1)) {
            const std::string prefix = path.substr(0, dot);
            const auto info = core::fieldInfo(project, prefix);
            if (info && info->kind == core::FieldKind::Records &&
                std::find(out.begin(), out.end(), prefix) == out.end()) {
                out.push_back(prefix);
            }
        }
        out.push_back(path);
    }
    std::vector<std::string> lists;
    recordLists(project, core::toPlain(project), {}, lists);
    for (const std::string& list : lists) {  // the empty ones
        if (std::find(out.begin(), out.end(), list) == out.end()) {
            out.push_back(list);
        }
    }
    return out;
}

}  // namespace mdbridge::app::forms
