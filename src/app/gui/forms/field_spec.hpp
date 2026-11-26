#pragma once
// How a field of the project is edited (milestone 5 W0-C): the kind of widget, its range hint,
// choices, file filters and texts. defaultFieldSpec() gives MDATK v0.0.4.1's choice for each field
// its GUI binds (mdatk/gui/app.py); pages may start from it and change labels or choices.

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <QString>
#include <QStringList>

#include "mdbridge/core/value.hpp"

namespace mdbridge::app::forms {

// The widget that edits a field (Auto: chosen from the field's kind, see createFieldWidget()).
enum class Editor {
    Auto,
    Number,         // spin box: Integer and Number fields
    Line,           // one-line text: text, optional numbers, comma-separated lists, JSON
    Path,           // line + Browse (file dialog)
    Combo,          // fixed choices; an unknown value is kept and shown (AP-05)
    EditableCombo,  // free text with suggestions
    Color,          // line + swatch + picker; any Matplotlib colour text is kept (AP-03)
    Check,          // check box (Boolean, OptionalBoolean)
    MultiLine,      // multi-line text
    Table,          // list editor (QTableView): lists of records, of text, of number rows
};

// What a file dialog picks.
enum class PathMode { OpenFile, OpenFiles, SaveFile, Folder };

struct Choice {
    QString value;  // the text stored
    QString label;  // the text shown ("" = value)
};

// A file dialog filter: MDATK's (name, patterns), e.g. ("Topology files", {"*.parm7", ...}).
struct FileFilter {
    QString name;
    QStringList patterns;
};

struct FieldSpec;

// A column of a list editor: a key of the record ("" for a list of plain items, "0", "1", ... for
// rows of numbers).
struct Column {
    std::string key;
    QString title;
};

struct FieldSpec {
    Editor editor = Editor::Auto;
    QString label;        // short name: accessible name, check box text ("" = from the field's location)
    QString hint;         // tooltip (shown below a problem's message)
    QString placeholder;  // shown while the text is empty ("last", "from file", "all groups")
    // Number: MDATK's spin range, shown as a hint; values outside are kept and flagged, never clamped
    std::optional<double> minimum;
    std::optional<double> maximum;
    double step = 1.0;
    int decimals = -1;  // stepped values are rounded to this many decimals (-1: integers / none)
    // Combo, EditableCombo (and list editor columns of either kind)
    std::vector<Choice> choices;
    // Path (and the Browse button of a list editor)
    PathMode pathMode = PathMode::OpenFile;
    QString dialogTitle;
    std::vector<FileFilter> filters;
    // Color: "" is a valid value (output.charts.color: the palette)
    bool allowEmpty = false;
    // MultiLine
    bool monospace = false;
    // Table: the columns (empty: from the record's fields), the item an Add inserts (Null: the
    // record's defaults, "" or zeros), the record key a browsed file goes to ("" for plain lists)
    std::vector<Column> columns;
    core::Value newItem;
    std::string browseKey;
    bool browse = false;  // a Browse… button that appends the files picked (pathMode, filters)
};

// MDATK v0.0.4.1's widget for the field at `modelPath` (the path without list indices,
// core::FieldInfo::modelPath): spin ranges, steps and digits, combos and their labels, file
// chooser filters and titles, placeholders, colours, multi-line texts, list editors. A field
// MDATK's GUI does not bind gets the default spec (Editor::Auto, no texts).
FieldSpec defaultFieldSpec(std::string_view modelPath);

// "Topology files (*.parm7 *.PARM7 *.prmtop ...)": the Qt name filter of a file filter. MDATK's
// filters match in any case; Qt's file dialogs match case-sensitively on Linux, so each pattern is
// also given in lower and upper case.
QString nameFilter(const FileFilter& filter);
// The name filters of a spec, followed by "All files (*)" when there is at least one (MDATK's
// add_filters).
QStringList nameFilters(const FieldSpec& spec);

// MDATK's filters (app.py).
FileFilter topologyFiles();
FileFilter trajectoryFiles();
FileFilter referenceFiles();
FileFilter dataFiles();
FileFilter mapFiles();
FileFilter projectFiles();
FileFilter runLogFiles();

}  // namespace mdbridge::app::forms
