#include "gui_harness.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>

#include <QAbstractButton>
#include <QAbstractSlider>
#include <QAbstractSpinBox>
#include <QApplication>
#include <QComboBox>
#include <QCoreApplication>
#include <QFont>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QMouseEvent>
#include <QPixmap>
#include <QPlainTextEdit>
#include <QStringList>
#include <QTableView>
#include <QTest>
#include <QTextEdit>
#include <QWheelEvent>
#include <QWidget>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/main_window.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_file.hpp"
#include "mdbridge/core/settings.hpp"
#include "mdbridge/render/text.hpp"
#include "mdbridge/render/text_qt.hpp"
#include "gui/pages/layout_page.hpp"
#include "gui/preview_widget.hpp"

namespace fs = std::filesystem;

namespace mdbridge::guitest {

namespace {

constexpr int kDefaultWaitMs = 60000;

QString toQString(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw CaseError("cannot read " + core::pathToUtf8(path));
    }
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string genericUtf8(const fs::path& path) {
    const std::u8string text = path.generic_u8string();
    return std::string(text.begin(), text.end());
}

void replaceAll(std::string& text, const std::string& from, const std::string& to) {
    if (from.empty()) {
        return;
    }
    for (std::size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
}

const std::string& stringEntry(const core::Value& table, std::string_view key, const std::string& where) {
    const core::Value* v = table.find(key);
    const std::string* s = v ? v->get<std::string>() : nullptr;
    if (!s) {
        throw CaseError(where + ": '" + std::string(key) + "' must be a string");
    }
    return *s;
}

std::string optionalString(const core::Value& table, std::string_view key) {
    const core::Value* v = table.find(key);
    return v ? v->stringOr() : std::string();
}

std::vector<std::pair<std::string, std::string>> fileMap(const core::Value& root, std::string_view key,
                                                         const std::string& where) {
    std::vector<std::pair<std::string, std::string>> out;
    const core::Value* v = root.find(key);
    if (!v) {
        return out;
    }
    const core::Table* table = v->get<core::Table>();
    if (!table) {
        throw CaseError(where + ": '" + std::string(key) + "' must be an object");
    }
    for (const auto& [name, file] : *table) {
        const std::string* text = file.get<std::string>();
        if (!text) {
            throw CaseError(where + ": '" + std::string(key) + "." + name + "' must be a path under tests/data");
        }
        out.emplace_back(name, *text);
    }
    return out;
}

int digits(const std::string& text, std::size_t at, std::size_t count, const std::string& where) {
    int value = 0;
    const char* first = text.data() + at;
    const auto [end, ec] = std::from_chars(first, first + count, value);
    if (ec != std::errc() || end != first + count) {
        throw CaseError(where + ": 'clock' must be YYYY-MM-DDTHH:MM:SS, not '" + text + "'");
    }
    return value;
}

core::DateTime parseClock(const std::string& text, const std::string& where) {
    if (text.size() != 19 || text[4] != '-' || text[7] != '-' || text[10] != 'T' || text[13] != ':' || text[16] != ':') {
        throw CaseError(where + ": 'clock' must be YYYY-MM-DDTHH:MM:SS, not '" + text + "'");
    }
    core::DateTime t;
    t.kind = core::DateTime::Kind::LocalDateTime;
    t.year = digits(text, 0, 4, where);
    t.month = digits(text, 5, 2, where);
    t.day = digits(text, 8, 2, where);
    t.hour = digits(text, 11, 2, where);
    t.minute = digits(text, 14, 2, where);
    t.second = digits(text, 17, 2, where);
    return t;
}

int daysInMonth(int year, int month) {
    static constexpr int kDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    const bool leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    return month == 2 && leap ? 29 : kDays[month - 1];
}

void addSecond(core::DateTime& t) {
    if (++t.second < 60) return;
    t.second = 0;
    if (++t.minute < 60) return;
    t.minute = 0;
    if (++t.hour < 24) return;
    t.hour = 0;
    if (++t.day <= daysInMonth(t.year, t.month)) return;
    t.day = 1;
    if (++t.month <= 12) return;
    t.month = 1;
    ++t.year;
}

// The text a user types for `value` (a list of strings: "a, b", as the comma-separated fields).
QString typedText(const core::Value& value, const char* widgetType) {
    if (const auto* s = value.get<std::string>()) return toQString(*s);
    if (const auto* items = value.get<core::Array>()) {
        QStringList parts;
        for (const core::Value& item : *items) {
            const auto* text = item.get<std::string>();
            if (!text) {
                throw CaseError(std::string("a ") + widgetType + " takes a list of strings only");
            }
            parts.append(toQString(*text));
        }
        return parts.join(QStringLiteral(", "));
    }
    if (const auto* i = value.get<std::int64_t>()) return QString::number(*i);
    if (const auto* d = value.get<double>()) return toQString(core::pythonFloatRepr(*d));
    throw CaseError(std::string("a ") + widgetType + " takes a string or a number, not " +
                    std::string(core::typeName(value.type())));
}

void typeInto(QWidget* widget, const QString& text) {
    if (text.isEmpty()) {
        QTest::keyClick(widget, Qt::Key_Delete);  // the selection: the field is cleared
        return;
    }
    const QStringList lines = text.split(QLatin1Char('\n'));
    for (qsizetype k = 0; k < lines.size(); ++k) {
        if (k > 0) QTest::keyClick(widget, Qt::Key_Return);
        QTest::keyClicks(widget, lines[k]);
    }
}

// Select all, type `text` ("" clears), Enter: line edits, spin boxes, editable combo boxes.
void typeLine(QWidget* widget, const QString& text) {
    QWidget* target = widget;
    if (auto* box = qobject_cast<QComboBox*>(widget)) {
        box->lineEdit()->selectAll();
        target = box->lineEdit();
    } else if (auto* spin = qobject_cast<QAbstractSpinBox*>(widget)) {
        spin->selectAll();
    } else if (auto* edit = qobject_cast<QLineEdit*>(widget)) {
        edit->selectAll();
    }
    typeInto(target, text);
    QTest::keyClick(target, Qt::Key_Return);
}

// The item of a drop-down whose project value (item data) is `text`.
void chooseItem(QComboBox* box, const QString& text, const std::string& name) {
    const int index = box->findData(text);
    if (index < 0) {
        throw CaseError(name + ": the drop-down has no item whose value is '" + text.toStdString() + "'");
    }
    box->setCurrentIndex(index);
    emit box->activated(index);  // as when a user picks the item
}

// Multi-line text is stored when the field is left: its field widget commits as then.
void commitAsOnLeaving(QWidget* widget) {
    for (QWidget* w = widget; w; w = w->parentWidget()) {
        if (auto* field = qobject_cast<app::forms::FieldWidget*>(w)) {
            field->commit();
            return;
        }
    }
    widget->clearFocus();
}

QList<QWidget*> widgetsNamed(QWidget* root, const QString& key) {
    QList<QWidget*> found = root->findChildren<QWidget*>(key);
    if (root->objectName() == key) {
        found.prepend(root);
    }
    return found;
}

void processEvents() {
    QApplication::processEvents();
    QApplication::processEvents();  // the view stores a closed editor's data queued
}

// Opens the cell's editor, `enter`s into it, presses Enter and waits until the editor is closed.
template <class Enter>
void editCell(app::forms::ListEditor* list, const std::string& key, Enter&& enter) {
    QTableView* view = list->view();
    const QModelIndex index = list->cellIndex(key);
    view->setFocus();
    view->setCurrentIndex(index);
    view->edit(index);
    QWidget* editor = view->indexWidget(index);
    if (!editor) {
        throw CaseError(key + ": the table cell opened no editor");
    }
    enter(editor);
    QTest::keyClick(editor, Qt::Key_Return);  // the delegate's filter on the editor commits and closes
    if (!QTest::qWaitFor([&] { return !view->indexWidget(list->cellIndex(key)); }, 5000)) {
        throw CaseError(key + ": the table cell's editor did not close on Enter");
    }
    processEvents();
}

void enterCellText(QWidget* editor, const QString& text, const std::string& key) {
    if (auto* box = qobject_cast<QComboBox*>(editor)) {
        if (!box->isEditable()) {
            const int index = box->findData(text);
            if (index < 0) {
                throw CaseError(key + ": the drop-down has no item whose value is '" + text.toStdString() + "'");
            }
            box->setCurrentIndex(index);
            return;
        }
        box->lineEdit()->selectAll();
        typeInto(box->lineEdit(), text);
        return;
    }
    if (auto* edit = qobject_cast<QLineEdit*>(editor)) {
        edit->selectAll();
        typeInto(edit, text);
        return;
    }
    throw CaseError(key + ": cannot type into a " + std::string(editor->metaObject()->className()));
}

// A check cell: toggled with Space when its state differs.
void toggleCell(app::forms::ListEditor* list, const std::string& key, bool wanted) {
    const auto checked = [&] { return list->cellIndex(key).data(Qt::CheckStateRole).toInt() == Qt::Checked; };
    if (checked() != wanted) {
        list->view()->setFocus();
        list->view()->setCurrentIndex(list->cellIndex(key));
        QTest::keyClick(list->view(), Qt::Key_Space);
        processEvents();
    }
    if (checked() != wanted) {
        throw CaseError(key + ": pressing Space did not change the check cell");
    }
}

// The cell of field path `key`; throws when neither a widget nor a cell has it.
app::forms::CellLocation cellOf(QWidget* root, const QString& key) {
    const app::forms::CellLocation cell = app::forms::findCell(root, key.toStdString());
    if (!cell.editor) {
        throw CaseError("no widget or table cell is named '" + key.toStdString() + "'");
    }
    return cell;
}

}  // namespace

std::vector<fs::path> findCases(const fs::path& dir) {
    std::vector<fs::path> out;
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) {
        return out;
    }
    for (const fs::directory_entry& e : fs::recursive_directory_iterator(dir)) {
        if (e.is_regular_file() && e.path().extension() == ".json") {
            out.push_back(e.path());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

GuiCase loadCase(const fs::path& file, const fs::path& casesDir) {
    GuiCase c;
    c.source = file;
    c.name = genericUtf8(file.lexically_relative(casesDir).replace_extension());
    const std::string where = c.name;
    core::Value root;
    try {
        root = core::parseJson(readFile(file), core::pathToUtf8(file));
    } catch (const std::exception& e) {
        throw CaseError(e.what());
    }
    if (!root.isTable()) {
        throw CaseError(where + ": a case is a JSON object");
    }
    if (optionalString(root, "format") != "mdbridge-gui-case") {
        throw CaseError(where + ": 'format' must be \"mdbridge-gui-case\"");
    }
    const core::Value* version = root.find("version");
    if (!version || !version->get<std::int64_t>() || *version->get<std::int64_t>() != 1) {
        throw CaseError(where + ": 'version' must be 1 (this harness reads version 1)");
    }
    c.description = optionalString(root, "description");
    if (root.find("clock")) {
        c.clock = parseClock(stringEntry(root, "clock", where), where);
    }
    c.files = fileMap(root, "files", where);
    c.expected = fileMap(root, "expected", where);
    const core::Value* actions = root.find("actions");
    if (!actions || !actions->isArray()) {
        throw CaseError(where + ": 'actions' must be an array");
    }
    if (const core::Value* needs = root.find("needs")) {
        if (!needs->isArray()) {
            throw CaseError(where + ": 'needs' must be an array of strings");
        }
        for (const core::Value& need : *needs->get<core::Array>()) {
            if (!need.isString()) {
                throw CaseError(where + ": 'needs' must be an array of strings");
            }
            c.needs.push_back(need.stringOr());
        }
    }
    if (const core::Value* start = root.find("start")) {
        if (!start->isTable() || !start->find("project") || !start->find("project")->isString()) {
            throw CaseError(where + ": 'start' must be {\"project\": name}");
        }
        c.start = start->find("project")->stringOr();
    }
    c.settings = optionalString(root, "settings");
    if (const core::Value* log = root.find("log")) {
        if (!log->isArray()) throw CaseError(where + ": 'log' must be an array of strings");
        std::vector<std::string> lines;
        for (const core::Value& line : *log->get<core::Array>()) {
            if (!line.isString()) throw CaseError(where + ": 'log' must be an array of strings");
            lines.push_back(line.stringOr());
        }
        c.log = std::move(lines);
    }
    if (const core::Value* outcome = root.find("outcome")) {
        if (!outcome->isTable()) throw CaseError(where + ": 'outcome' must be an object");
        c.outcome = *outcome;
    }
    if (const core::Value* canvas = root.find("canvas")) {
        const core::Value* w = canvas->isTable() ? canvas->find("width") : nullptr;
        const core::Value* h = canvas->isTable() ? canvas->find("height") : nullptr;
        if (!w || !h || !w->get<std::int64_t>() || !h->get<std::int64_t>()) {
            throw CaseError(where + ": 'canvas' must be {\"width\": W, \"height\": H}");
        }
        c.canvas = std::pair<int, int>(static_cast<int>(*w->get<std::int64_t>()), static_cast<int>(*h->get<std::int64_t>()));
    }
    c.actions = *actions->get<core::Array>();
    for (std::size_t k = 0; k < c.actions.size(); ++k) {
        if (!c.actions[k].isTable() || !c.actions[k].find("action") || !c.actions[k].find("action")->isString()) {
            throw CaseError(where + ": action " + std::to_string(k + 1) + " must be an object with an \"action\" string");
        }
    }
    return c;
}

bool GuiCase::needsAny(std::string_view need) const {
    return std::find(needs.begin(), needs.end(), need) != needs.end();
}

std::vector<std::string> unsupportedNeeds(const GuiCase& c) {
    static const std::vector<std::string> supported{"answers", "canvas", "close", "lossless", "new", "open_failure",
                                                    "save_as", "settings", "snapshot", "start", "table",
                                                    "toggle", "type"};
    std::vector<std::string> out;
    for (const std::string& need : c.needs) {
        if (std::find(supported.begin(), supported.end(), need) == supported.end()) {
            out.push_back(need);
        }
    }
    return out;
}

FixedClock::FixedClock(const core::DateTime& start) {
    auto next = std::make_shared<core::DateTime>(start);
    m_previous = core::setClock(core::Clock{[next] {
                                                const core::DateTime now = *next;
                                                addSecond(*next);
                                                return now;
                                            },
                                            {}});
}

FixedClock::~FixedClock() { core::setClock(m_previous); }

QWidget* findField(QWidget* root, const QString& key) {
    QList<QWidget*> found = root->findChildren<QWidget*>(key);
    if (root->objectName() == key) {
        found.prepend(root);
    }
    if (found.isEmpty()) {
        throw CaseError("no widget is named '" + key.toStdString() + "'");
    }
    if (found.size() > 1) {
        throw CaseError(std::to_string(found.size()) + " widgets are named '" + key.toStdString() + "'");
    }
    return found.front();
}

void setField(QWidget* widget, const core::Value& value) {
    const std::string name = widget->objectName().toStdString();
    if (auto* box = qobject_cast<QComboBox*>(widget)) {
        const QString text = typedText(value, "combo box");
        if (box->isEditable()) {
            box->lineEdit()->selectAll();
            QTest::keyClicks(box->lineEdit(), text);
            QTest::keyClick(box->lineEdit(), Qt::Key_Return);
            return;
        }
        chooseItem(box, text, name);
        return;
    }
    if (auto* edit = qobject_cast<QLineEdit*>(widget)) {
        edit->selectAll();
        QTest::keyClicks(edit, typedText(value, "line edit"));
        QTest::keyClick(edit, Qt::Key_Return);
        return;
    }
    if (auto* spin = qobject_cast<QAbstractSpinBox*>(widget)) {
        spin->selectAll();
        QTest::keyClicks(spin, typedText(value, "spin box"));
        QTest::keyClick(spin, Qt::Key_Return);
        return;
    }
    if (auto* plain = qobject_cast<QPlainTextEdit*>(widget)) {
        plain->selectAll();
        typeInto(plain, typedText(value, "text edit"));
        return;
    }
    if (auto* rich = qobject_cast<QTextEdit*>(widget)) {
        rich->selectAll();
        typeInto(rich, typedText(value, "text edit"));
        return;
    }
    if (auto* button = qobject_cast<QAbstractButton*>(widget); button && button->isCheckable()) {
        const bool* wanted = value.get<bool>();
        if (!wanted) {
            throw CaseError(name + ": a check box takes true or false");
        }
        if (button->isChecked() != *wanted) {
            QTest::keyClick(button, Qt::Key_Space);  // a click needs the indicator's place, which depends on the style
        }
        if (button->isChecked() != *wanted) {
            throw CaseError(name + ": pressing Space did not change the check box");
        }
        return;
    }
    if (auto* slider = qobject_cast<QAbstractSlider*>(widget)) {
        const auto* d = value.get<double>();
        const auto* i = value.get<std::int64_t>();
        if (!d && !i) {
            throw CaseError(name + ": a slider takes a number");
        }
        const double scale = slider->property("mdbridgeScale").isValid() ? slider->property("mdbridgeScale").toDouble() : 1.0;
        const double v = (d ? *d : static_cast<double>(*i)) * scale;
        const int units = static_cast<int>(std::lround(v));
        if (units < slider->minimum() || units > slider->maximum()) {
            throw CaseError(name + ": the value is outside the slider's range");
        }
        slider->setFocus();
        slider->setValue(units);  // as a drag to there
        return;
    }
    if (const auto* items = value.get<core::Array>()) {
        QList<QAbstractButton*> choices;
        for (QAbstractButton* b : widget->findChildren<QAbstractButton*>()) {
            if (b->property("mdbridgeChoice").isValid()) {
                choices.append(b);
            }
        }
        if (!choices.isEmpty()) {
            QStringList wanted;
            for (const core::Value& item : *items) {
                wanted.append(typedText(item, "check box list"));
            }
            for (QAbstractButton* b : choices) {
                if (b->isChecked() && !wanted.contains(b->property("mdbridgeChoice").toString())) {
                    QTest::keyClick(b, Qt::Key_Space);
                }
            }
            for (const QString& w : wanted) {
                QAbstractButton* box = nullptr;
                for (QAbstractButton* b : choices) {
                    if (b->property("mdbridgeChoice").toString() == w) box = b;
                }
                if (!box) {
                    throw CaseError(name + ": there is no check box for '" + w.toStdString() + "'");
                }
                if (!box->isChecked()) {
                    QTest::keyClick(box, Qt::Key_Space);
                }
            }
            return;
        }
        for (std::size_t k = 0; k < items->size(); ++k) {
            const QString part = widget->objectName() + QStringLiteral(".%1").arg(k);
            const QList<QWidget*> found = widget->findChildren<QWidget*>(part);
            if (found.size() != 1) {
                throw CaseError(name + ": no single part '" + part.toStdString() + "' takes item " + std::to_string(k));
            }
            setField(found.front(), (*items)[k]);
            QApplication::processEvents();
        }
        return;
    }
    throw CaseError(name + ": cannot enter a value into a " + std::string(widget->metaObject()->className()));
}

void typeField(QWidget* widget, const QString& text) {
    const std::string name = widget->objectName().toStdString();
    if (qobject_cast<QLineEdit*>(widget) || qobject_cast<QAbstractSpinBox*>(widget)) {
        typeLine(widget, text);
        return;
    }
    if (auto* box = qobject_cast<QComboBox*>(widget); box && box->isEditable()) {
        typeLine(widget, text);
        return;
    }
    if (auto* plain = qobject_cast<QPlainTextEdit*>(widget)) {
        plain->setFocus();
        plain->selectAll();
        typeInto(plain, text);
        commitAsOnLeaving(plain);
        return;
    }
    if (auto* rich = qobject_cast<QTextEdit*>(widget)) {
        rich->setFocus();
        rich->selectAll();
        typeInto(rich, text);
        commitAsOnLeaving(rich);
        return;
    }
    throw CaseError(name + ": cannot type into a " + std::string(widget->metaObject()->className()));
}

void setKey(QWidget* root, const QString& key, const core::Value& value) {
    if (!widgetsNamed(root, key).isEmpty()) {
        setField(findField(root, key), value);
        processEvents();
        return;
    }
    const app::forms::CellLocation cell = cellOf(root, key);
    const std::string name = key.toStdString();
    if (cell.index.data(Qt::CheckStateRole).isValid()) {
        const bool* wanted = value.get<bool>();
        if (!wanted) {
            throw CaseError(name + ": a check cell takes true or false");
        }
        toggleCell(cell.editor, name, *wanted);
        return;
    }
    const QString text = typedText(value, "table cell");
    editCell(cell.editor, name, [&](QWidget* editor) { enterCellText(editor, text, name); });
}

void typeKey(QWidget* root, const QString& key, const QString& text) {
    if (!widgetsNamed(root, key).isEmpty()) {
        typeField(findField(root, key), text);
        processEvents();
        return;
    }
    const app::forms::CellLocation cell = cellOf(root, key);
    const std::string name = key.toStdString();
    if (cell.index.data(Qt::CheckStateRole).isValid()) {
        throw CaseError(name + ": a check cell takes no text (use \"set\")");
    }
    editCell(cell.editor, name, [&](QWidget* editor) { enterCellText(editor, text, name); });
}

void tableAction(QWidget* root, const QString& key, const std::string& what, int row, int by) {
    app::forms::ListEditor* list = nullptr;
    for (app::forms::ListEditor* e : root->findChildren<app::forms::ListEditor*>()) {
        if (e->path() == key.toStdString()) {
            list = e;
        }
    }
    const std::string name = key.toStdString();
    if (!list) {
        throw CaseError("no list editor edits '" + name + "'");
    }
    const auto click = [&](const char* button) {
        auto* b = qobject_cast<QAbstractButton*>(findField(root, key + QLatin1String(button)));
        if (!b) {
            throw CaseError(name + button + " is no button");
        }
        QTest::mouseClick(b, Qt::LeftButton);
        processEvents();
    };
    const auto select = [&] {
        if (row < 0 || row >= list->rowCount()) {
            throw CaseError(name + ": there is no row " + std::to_string(row) + " (" + std::to_string(list->rowCount()) +
                            " rows)");
        }
        list->selectRow(row);
    };
    if (what == "add") {
        click(":add");
    } else if (what == "remove") {
        select();
        click(":remove");
    } else if (what == "move") {
        if (by == 0) {
            throw CaseError(name + ": 'by' must be -1 or +1");
        }
        select();
        for (int k = 0; k < (by < 0 ? -by : by); ++k) {
            click(by < 0 ? ":up" : ":down");
        }
    } else {
        throw CaseError(name + ": 'do' must be \"add\", \"remove\" or \"move\", not '" + what + "'");
    }
}

fs::path saveScreenshot(QWidget* widget, const fs::path& dir, const std::string& name) {
    std::error_code ec;
    fs::create_directories(dir, ec);
    const fs::path file = dir / core::pathFromUtf8(name + ".png");
    const QPixmap pixmap = widget->grab();
    if (pixmap.isNull() || !pixmap.save(toQString(core::pathToUtf8(file)), "PNG")) {
        throw CaseError("cannot save the screenshot " + core::pathToUtf8(file));
    }
    return file;
}

void useBundledFont() {
    const std::string_view name = render::bundledFont(render::FontFace::Regular).family;
    const QString family = QString::fromUtf8(name.data(), static_cast<qsizetype>(name.size()));
    QStringList families;
    try {
        families = render::registerBundledFonts();
    } catch (const std::exception& e) {
        qWarning("the bundled fonts: %s", e.what());  // the layout tests then report what they measure
    }
    if (!families.contains(family)) {
        qWarning("the bundled fonts have no %s: the tests run with Qt's default font", qPrintable(family));
        return;
    }
    QApplication::setFont(render::interfaceFont(9));  // as the application sets it, at 9 pt
}

fs::path screenshotDir(const fs::path& fallback) {
    if (const char* dir = std::getenv("MDBRIDGE_GUI_SCREENSHOTS"); dir && *dir) {
        return core::pathFromUtf8(dir);
    }
    return fallback;
}

std::string replaceFolder(std::string text, const std::string& folder, const std::string& name, bool anySeparator) {
    if (folder.empty()) {
        return text;
    }
    if (!anySeparator) {
        replaceAll(text, folder, name);
        return text;
    }
    const auto slash = [](char c) { return c == '\\' ? '/' : c; };
    // where a path in a message ends: "Could not open <path>: <reason>", "<path> (save them ...)"
    const auto endsPath = [](char c) {
        return c == ' ' || c == '\t' || c == '\n' || c == '"' || c == '\'' || c == ':' || c == ',' || c == ';' ||
               c == ')' || c == ']' || c == '}';
    };
    std::string out;
    out.reserve(text.size());
    std::size_t at = 0;
    while (at < text.size()) {
        bool found = at + folder.size() <= text.size();
        for (std::size_t k = 0; found && k < folder.size(); ++k) {
            found = slash(text[at + k]) == slash(folder[k]);
        }
        if (!found) {
            out += text[at++];
            continue;
        }
        out += name;
        for (at += folder.size(); at < text.size() && !endsPath(text[at]); ++at) {
            out += slash(text[at]);
        }
    }
    return out;
}

namespace {

// Windows: the harness's folders are Qt's '/' paths joined by std::filesystem ('\\'); the window
// shows a path as it was given (mixed) or normalised by core::ntpath ('/'), and the cases were
// recorded with '/': the folders match whatever their separators and the paths after them
// compare with '/'.
#ifdef _WIN32
constexpr bool kAnySeparator = true;
#else
constexpr bool kAnySeparator = false;
#endif

std::string withFolder(std::string text, const fs::path& dir, const std::string& name) {
    const fs::path normal = dir.lexically_normal();
    text = replaceFolder(std::move(text), genericUtf8(normal), name, kAnySeparator);
    return replaceFolder(std::move(text), core::pathToUtf8(normal), name, kAnySeparator);
}

}  // namespace

std::string normalise(std::string text, const fs::path& caseDir) {
    replaceAll(text, "\r\n", "\n");
    if (!caseDir.empty()) {
        text = withFolder(std::move(text), caseDir, "{case_dir}");
    }
    return text;
}

std::string compareText(const std::string& actual, const std::string& expected) {
    if (actual == expected) {
        return {};
    }
    std::istringstream a(actual);
    std::istringstream b(expected);
    std::string la;
    std::string lb;
    for (int line = 1;; ++line) {
        const bool hasA = static_cast<bool>(std::getline(a, la));
        const bool hasB = static_cast<bool>(std::getline(b, lb));
        if (!hasA && !hasB) {
            return "the texts differ in their final line end";
        }
        if (!hasA || !hasB || la != lb) {
            std::ostringstream out;
            out << "line " << line << ": expected " << (hasB ? "'" + lb + "'" : std::string("the end of the file"))
                << ", got " << (hasA ? "'" + la + "'" : std::string("the end of the file"));
            return out.str();
        }
    }
}

CaseRunner::CaseRunner(app::MainWindow& window, fs::path caseDir, fs::path screenshots, std::string caseName)
    : m_window(window), m_caseDir(std::move(caseDir)), m_screenshotDir(std::move(screenshots)),
      m_caseName(std::move(caseName)) {
    std::replace(m_caseName.begin(), m_caseName.end(), '/', '_');
    m_window.setPrompter(&m_prompter);
    m_previousDialog = app::forms::setFileDialogFunction(
        [this](QWidget*, const app::forms::FileDialogRequest& request) { return dialogAnswer(request); });
}

CaseRunner::~CaseRunner() {
    app::forms::setFileDialogFunction(m_previousDialog);
    m_window.setPrompter(nullptr);
}

// The file dialogs of an action: Save as… gets `name`, Open `path`, Save of a project without a
// file answers.name (all relative to the case folder); any other dialog is cancelled.
QStringList CaseRunner::dialogAnswer(const app::forms::FileDialogRequest& request) const {
    const std::string what = m_current.isTable() ? optionalString(m_current, "action") : std::string();
    std::string name;
    if (what == "save_as") name = optionalString(m_current, "name");
    else if (what == "open") name = optionalString(m_current, "path");
    else if (const core::Value* answers = m_current.isTable() ? m_current.find("answers") : nullptr)
        name = optionalString(*answers, "name");
    else if (what == "save") name = optionalString(m_current, "path");
    (void)request;
    if (name.empty()) return {};
    return {toQString(core::pathToUtf8(m_caseDir / core::pathFromUtf8(name)))};
}

void CaseRunner::start(const std::string& project) {
    m_index = -1;
    m_prompter.action = -1;
    m_prompter.answers = core::Value();
    m_window.startWith(m_caseDir / core::pathFromUtf8(project));
    waitForPreview();
}

// As the generator's flush() after the window is shown: the first drawing is done (and its
// warnings known) before the first action.
void CaseRunner::waitForPreview() {
    if (!QTest::qWaitFor([this] { return m_window.isHidden() || m_window.previewSettled(); }, 30000)) {
        throw CaseError("the preview did not settle in 30 s");
    }
}

void CaseRunner::run(const core::Array& actions) {
    waitForPreview();
    for (std::size_t k = 0; k < actions.size(); ++k) {
        m_index = static_cast<int>(k);
        m_current = actions[k];
        m_prompter.action = m_index;
        const core::Value* answers = actions[k].find("answers");
        m_prompter.answers = answers ? *answers : core::Value();
        try {
            runAction(actions[k]);
            // as the generator's flush(): the redraws the action caused are done (and their warnings
            // logged) before the next one
            if (!QTest::qWaitFor([this] { return m_window.isHidden() || m_window.previewSettled(); }, 30000)) {
                throw CaseError("the preview did not settle in 30 s");
            }
        } catch (const std::exception& e) {
            throw CaseError("action " + std::to_string(k + 1) + " (" + actions[k].find("action")->stringOr() +
                            "): " + e.what());
        }
    }
}

void CaseRunner::waitUntilLoaded(int timeoutMs) {
    const bool done = QTest::qWaitFor([this] { return !m_window.loading() && m_window.preview()->settled(); },
                                      timeoutMs);
    if (!done) {
        throw CaseError("the project was not loaded in " + std::to_string(timeoutMs / 1000) + " s");
    }
}

void CaseRunner::runAction(const core::Value& action) {
    const std::string& what = stringEntry(action, "action", "action");
    const auto showPageOf = [this](const std::string& key) { m_window.revealField(toQString(key)); };
    if (what == "open") {
        const std::size_t dialogs = m_prompter.dialogs.size();
        if (m_current.isNull()) m_current = action;  // runAction called on its own
        m_window.chooseAndOpen();
        if (m_prompter.dialogs.size() > dialogs && !m_openMayFail) {
            throw CaseError("cannot open " + stringEntry(action, "path", what) + ": " +
                            m_prompter.dialogs.back().find("detail")->stringOr());
        }
    } else if (what == "new") {
        m_window.newProject();
    } else if (what == "close") {
        m_window.close();
    } else if (what == "save_as") {
        if (m_current.isNull()) m_current = action;
        m_window.saveAs();
    } else if (what == "snapshot") {
        snapshot(action);
    } else if (what == "wait") {
        int timeoutMs = kDefaultWaitMs;
        if (const core::Value* t = action.find("timeout_s")) {
            if (const auto* i = t->get<std::int64_t>()) timeoutMs = static_cast<int>(*i * 1000);
            else if (const auto* d = t->get<double>()) timeoutMs = static_cast<int>(*d * 1000.0);
            else throw CaseError("'timeout_s' must be a number");
        }
        waitUntilLoaded(timeoutMs);
    } else if (what == "set" || what == "toggle") {
        const core::Value* value = action.find("value");
        if (!value) {
            throw CaseError("'value' is missing");
        }
        showPageOf(stringEntry(action, "key", what));
        setKey(&m_window, toQString(stringEntry(action, "key", what)), *value);
    } else if (what == "type") {
        showPageOf(stringEntry(action, "key", what));
        typeKey(&m_window, toQString(stringEntry(action, "key", what)), toQString(stringEntry(action, "text", what)));
    } else if (what == "table") {
        const auto integer = [&](std::string_view k, int fallback) {
            const core::Value* v = action.find(k);
            if (!v) return fallback;
            const auto* i = v->get<std::int64_t>();
            if (!i) throw CaseError("'" + std::string(k) + "' must be a whole number");
            return static_cast<int>(*i);
        };
        showPageOf(stringEntry(action, "key", what));
        tableAction(&m_window, toQString(stringEntry(action, "key", what)), stringEntry(action, "do", what),
                    integer("row", -1), integer("by", 0));
    } else if (what == "click") {
        QWidget* widget = findField(&m_window, toQString(stringEntry(action, "key", what)));
        auto* button = qobject_cast<QAbstractButton*>(widget);
        if (!button) {
            throw CaseError(widget->objectName().toStdString() + " is no button");
        }
        m_window.revealField(widget->objectName());
        QTest::mouseClick(button, Qt::LeftButton);
    } else if (what == "save") {
        // Save (Ctrl+S); with "path" (the harness's own cases) Save as… onto that file
        if (m_current.isNull()) m_current = action;
        if (optionalString(action, "path").empty()) m_window.save();
        else m_window.saveAs();
    } else if (what.rfind("canvas.", 0) == 0) {
        canvasAction(action, what);
    } else if (what == "screenshot") {
        const std::string widgetName = optionalString(action, "widget");
        QWidget* widget = widgetName.empty() ? static_cast<QWidget*>(&m_window) : findField(&m_window, toQString(widgetName));
        m_shots.push_back(saveScreenshot(widget, m_screenshotDir, m_caseName + "-" + stringEntry(action, "name", what)));
    } else {
        throw CaseError("unknown action (tests/data/README.md, \"gui\", lists the actions)");
    }
}

// ---- canvas actions ---------------------------------------------------------------------------

namespace {

double number(const core::Value& v) {
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    if (const auto* d = v.get<double>()) return *d;
    throw CaseError("a number was expected");
}

void sendMouse(QWidget* widget, QEvent::Type type, const QPointF& at, Qt::MouseButton button, Qt::MouseButtons buttons,
               Qt::KeyboardModifiers modifiers) {
    QMouseEvent event(type, at, widget->mapToGlobal(at), button, buttons, modifiers);
    QCoreApplication::sendEvent(widget, &event);
}

}  // namespace

void CaseRunner::canvasAction(const core::Value& action, const std::string& what) {
    app::PreviewWidget* canvas = m_window.preview();
    if (!QTest::qWaitFor([&] { return m_window.previewSettled() && canvas->current(); }, 30000)) {
        throw CaseError("the canvas shows no figure");
    }
    const std::vector<QRectF> rects = canvas->slotRects();
    // MDATK's record of this action, if the case has one
    const core::Value* record = nullptr;
    if (const auto* records = m_canvasRecords.get<core::Array>()) {
        for (const core::Value& r : *records) {
            const core::Value* a = r.find("action");
            if (a && a->get<std::int64_t>() && *a->get<std::int64_t>() == m_index) record = &r;
        }
    }
    if (record) {
        const core::Value* hit = record->find("slots");
        const auto* list = hit ? hit->get<core::Array>() : nullptr;
        if (!list || list->size() != rects.size()) {
            throw CaseError("the canvas has " + std::to_string(rects.size()) + " slots, MDATK's had " +
                            std::to_string(list ? list->size() : 0));
        }
        const core::Value* structure = record->find("structure_slot");
        if (structure && structure->get<std::int64_t>() && *structure->get<std::int64_t>() != canvas->structureSlot()) {
            throw CaseError("the structure panel is in slot " + std::to_string(canvas->structureSlot()) + ", MDATK's in " +
                            std::to_string(*structure->get<std::int64_t>()));
        }
        for (std::size_t k = 0; k < rects.size(); ++k) {
            const core::Array& r = *(*list)[k].get<core::Array>();
            const QRectF m(number(r[0]), number(r[1]), number(r[2]), number(r[3]));
            const QRectF& o = rects[k];
            const double d = std::max({std::abs(m.left() - o.left()), std::abs(m.top() - o.top()),
                                       std::abs(m.right() - o.right()), std::abs(m.bottom() - o.bottom())});
            if (d > 3.0) {
                std::ostringstream note;
                note << m_caseName << " action " << m_index << " slot " << k << ": hit area differs from MDATK's by " << d
                     << " px (MDATK " << m.x() << ", " << m.y() << ", " << m.width() << ", " << m.height() << "; here " << o.x()
                     << ", " << o.y() << ", " << o.width() << ", " << o.height() << ")";
                m_geometryNotes.push_back(note.str());
            }
        }
    }
    const auto point = [&](const core::Value& spec) {
        const core::Value* slot = spec.find("slot");
        const core::Value* at = spec.find("at");
        if (!slot || !slot->get<std::int64_t>() || !at || !at->isArray() || at->get<core::Array>()->size() != 2) {
            throw CaseError("a point is {\"slot\": N, \"at\": [fx, fy]}");
        }
        const auto k = static_cast<std::size_t>(*slot->get<std::int64_t>());
        if (k >= rects.size()) throw CaseError("the canvas has no slot " + std::to_string(k));
        const QRectF& r = rects[k];
        const core::Array& f = *at->get<core::Array>();
        return QPointF(r.x() + number(f[0]) * r.width(), r.y() + number(f[1]) * r.height());
    };
    // the events: MDATK's positions, else resolved here as make_gui_cases.py resolves them
    struct Event {
        std::string kind;
        QPointF at;
        int notches = 0;
    };
    std::vector<Event> events;
    if (record) {
        for (const core::Value& e : *record->find("events")->get<core::Array>()) {
            const core::Array& a = *e.get<core::Array>();
            Event ev{a[0].stringOr(), QPointF(number(a[1]), number(a[2])), 0};
            if (a.size() > 3) ev.notches = static_cast<int>(number(a[3]));
            events.push_back(ev);
        }
    } else if (what == "canvas.wheel") {
        const QPointF at = point(*action.find("at"));
        for (const core::Value& n : *action.find("notches")->get<core::Array>()) {
            events.push_back({"scroll", at, static_cast<int>(number(n))});
        }
    } else if (what == "canvas.double_click") {
        const QPointF at = point(*action.find("at"));
        for (const char* kind : {"press", "release", "press", "double", "release"}) events.push_back({kind, at, 0});
    } else if (what == "canvas.drag") {
        const QPointF from = point(*action.find("from"));
        events.push_back({"press", from, 0});
        if (const core::Value* to = action.find("to")) {
            events.push_back({"motion", point(*to), 0});
        } else {
            for (const core::Value& d : *action.find("path")->get<core::Array>()) {
                const core::Array& xy = *d.get<core::Array>();
                events.push_back({"motion", from + QPointF(number(xy[0]), number(xy[1])), 0});
            }
        }
        events.push_back({"release", events.back().at, 0});
    } else {
        throw CaseError("unknown canvas action");
    }
    const core::Value* shift = action.find("shift");
    const Qt::KeyboardModifiers mods = shift && shift->get<bool>() && *shift->get<bool>() ? Qt::ShiftModifier : Qt::NoModifier;
    bool shotTaken = false;
    for (const Event& e : events) {
        if (e.kind == "press") {
            sendMouse(canvas, QEvent::MouseButtonPress, e.at, Qt::LeftButton, Qt::LeftButton, mods);
        } else if (e.kind == "double") {
            sendMouse(canvas, QEvent::MouseButtonDblClick, e.at, Qt::LeftButton, Qt::LeftButton, mods);
        } else if (e.kind == "motion") {
            sendMouse(canvas, QEvent::MouseMove, e.at, Qt::NoButton, Qt::LeftButton, mods);
        } else if (e.kind == "release") {
            if (!shotTaken && !canvas->dragMode().isEmpty()) {
                // a mid-drag frame (CI artifact): the highlight and label of a panel drag
                QCoreApplication::processEvents();
                m_shots.push_back(saveScreenshot(canvas, m_screenshotDir,
                                                 m_caseName + "-drag-" + std::to_string(m_index)));
                shotTaken = true;
            }
            sendMouse(canvas, QEvent::MouseButtonRelease, e.at, Qt::LeftButton, Qt::NoButton, mods);
        } else if (e.kind == "scroll") {
            QWheelEvent wheel(e.at, canvas->mapToGlobal(e.at), QPoint(), QPoint(0, 120 * e.notches), Qt::NoButton, mods,
                              Qt::NoScrollPhase, false);
            QCoreApplication::sendEvent(canvas, &wheel);
            if (!QTest::qWaitFor([this] { return m_window.previewSettled(); }, 30000)) {
                throw CaseError("the canvas did not settle after a wheel turn");
            }
        } else {
            throw CaseError("unknown canvas event '" + e.kind + "'");
        }
    }
}

void prepareCase(const GuiCase& c, const fs::path& dataDir, const fs::path& caseDir) {
    for (const auto& [name, source] : c.files) {
        const fs::path target = caseDir / core::pathFromUtf8(name);
        std::error_code ec;
        fs::create_directories(target.parent_path(), ec);
        fs::copy_file(dataDir / core::pathFromUtf8(source), target, fs::copy_options::overwrite_existing, ec);
        if (ec) {
            throw CaseError(c.name + ": cannot copy " + source + ": " + ec.message());
        }
    }
}

std::string checkExpected(const GuiCase& c, const fs::path& dataDir, const fs::path& caseDir) {
    std::string problems;
    for (const auto& [name, source] : c.expected) {
        const fs::path written = caseDir / core::pathFromUtf8(name);
        std::string problem;
        if (!fs::exists(written)) {
            problem = "was not written";
        } else {
            problem = compareText(normalise(readFile(written), caseDir), normalise(readFile(dataDir / core::pathFromUtf8(source)), {}));
        }
        if (!problem.empty()) {
            problems += name + " (expected " + source + "): " + problem + "\n";
        }
    }
    return problems;
}


// ---- recording -------------------------------------------------------------------------------

QString RecordingPrompter::ask(QWidget*, const app::Question& question) {
    QString answer = question.choices.empty() ? QString() : question.choices.front().first;
    if (answers.isTable()) {
        if (const core::Value* given = answers.find(question.kind.toStdString())) {
            answer = toQString(given->stringOr());
        }
    }
    bool known = false;
    for (const auto& [id, label] : question.choices) known = known || id == answer;
    if (!known) {
        throw CaseError("the answer '" + answer.toStdString() + "' is no choice of the " + question.kind.toStdString() +
                        " question");
    }
    questions.push_back(core::Value(core::Table{{"action", static_cast<std::int64_t>(action)},
                                                {"kind", question.kind.toStdString()},
                                                {"path", question.path.toStdString()},
                                                {"answer", answer.toStdString()}}));
    return answer;
}

void RecordingPrompter::message(QWidget*, app::MessageKind kind, const QString& text, const QString& detail) {
    dialogs.push_back(core::Value(core::Table{{"action", static_cast<std::int64_t>(action)},
                                              {"kind", kind == app::MessageKind::Error ? "error" : "info"},
                                              {"text", text.toStdString()},
                                              {"detail", detail.toStdString()}}));
}

void RecordingPrompter::printError(const QString& line) { stderrLines.push_back(line.toStdString()); }

void CaseRunner::snapshot(const core::Value& action) {
    core::Table rec{{"action", static_cast<std::int64_t>(m_index)}};
    app::forms::FormContext& context = m_window.formContext();
    context.commitPendingEdits();
    if (const core::Value* keys = action.find("keys"); keys && keys->isArray() && !keys->get<core::Array>()->empty()) {
        if (const QString errors = m_window.formErrorText(); !errors.isEmpty()) {
            rec.emplace_back("error", errors.toStdString());
        } else {
            core::Table values;
            for (const core::Value& key : *keys->get<core::Array>()) {
                values.emplace_back(key.stringOr(), m_window.document().get(key.stringOr()));
            }
            rec.emplace_back("values", core::Value(std::move(values)));
        }
    }
    if (const core::Value* editors = action.find("slots"); editors && editors->get<bool>() && *editors->get<bool>()) {
        app::pages::LayoutPage* layout = nullptr;
        for (const app::MainWindow::Page& page : m_window.pages()) {
            if (auto* l = qobject_cast<app::pages::LayoutPage*>(page.widget)) layout = l;
        }
        if (!layout) throw CaseError("the window has no Layout & panels page");
        core::Array out;
        const auto text = [](const QString& t) { return core::Value(t.toStdString()); };
        for (int k = 0; k < layout->slotCount(); ++k) {
            const app::pages::LayoutPage::SlotState st = layout->slotState(k);
            core::Array choices;
            for (const auto& [id, label] : st.choices) choices.emplace_back(core::Array{text(id), text(label)});
            out.emplace_back(core::Table{{"frame", text(st.frame)},
                                         {"type", text(st.type)},
                                         {"data", st.data ? text(*st.data) : core::Value()},
                                         {"choices", core::Value(std::move(choices))},
                                         {"group_a", text(st.groupA)},
                                         {"group_b", text(st.groupB)},
                                         {"group_a_enabled", st.groupAEnabled},
                                         {"group_b_enabled", st.groupBEnabled},
                                         {"more", text(st.more)},
                                         {"title", text(st.title)},
                                         {"options", text(st.options)}});
        }
        rec.emplace_back("slots", core::Value(std::move(out)));
    }
    if (const core::Value* page = action.find("page"); page && page->get<bool>() && *page->get<bool>()) {
        rec.emplace_back("page", m_window.currentPage().toStdString());
    }
    if (const core::Value* status = action.find("density_status"); status && status->get<bool>() && *status->get<bool>()) {
        auto* label = m_window.findChild<QLabel*>(QStringLiteral("density.status"));
        rec.emplace_back("density_status", label ? label->text().toStdString() : std::string("<no density.status label>"));
    }
    if (const core::Value* modified = action.find("modified"); modified && modified->get<bool>() && *modified->get<bool>()) {
        rec.emplace_back("modified", m_window.document().isModified());
    }
    m_snapshots.push_back(core::Value(std::move(rec)));
}

// ---- outcome ---------------------------------------------------------------------------------

std::string maskParserWords(std::string text) {
    // the TOML parser's own words differ from Python's tomllib (tests/core/test_project_file.cpp);
    // the file, line and column are compared
    static const std::string head = "not valid TOML: ";
    for (std::size_t at = text.find(head); at != std::string::npos; at = text.find(head, at + head.size())) {
        const std::size_t end = text.find(" (at line ", at + head.size());
        if (end == std::string::npos) break;
        text.replace(at + head.size(), end - at - head.size(), "<parser message>");
    }
    return text;
}

std::string placeholders(std::string text, const fs::path& caseDir, const fs::path& configDir) {
    replaceAll(text, "\r\n", "\n");
    text = maskParserWords(std::move(text));
    for (const auto& [dir, name] : {std::pair{configDir, "{config_dir}"}, std::pair{caseDir, "{case_dir}"}}) {
        if (dir.empty()) continue;
        text = withFolder(std::move(text), dir, name);
    }
    return text;
}

namespace {

double numberOf(const core::Value& v) {
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    return *v.get<double>();
}

bool isNumber(const core::Value& v) { return v.type() == core::Value::Type::Int || v.type() == core::Value::Type::Float; }

// a value with its texts' folders as placeholders
core::Value masked(const core::Value& v) {
    if (const auto* s = v.get<std::string>()) return core::Value(maskParserWords(*s));
    if (const auto* a = v.get<core::Array>()) {
        core::Array out;
        for (const core::Value& item : *a) out.push_back(masked(item));
        return core::Value(std::move(out));
    }
    if (const auto* t = v.get<core::Table>()) {
        core::Table out;
        for (const auto& [k, item] : *t) out.emplace_back(k, masked(item));
        return core::Value(std::move(out));
    }
    return v;
}

core::Value withPlaceholders(const core::Value& v, const fs::path& caseDir, const fs::path& configDir) {
    if (const auto* s = v.get<std::string>()) return core::Value(placeholders(*s, caseDir, configDir));
    if (const auto* a = v.get<core::Array>()) {
        core::Array out;
        for (const core::Value& item : *a) out.push_back(withPlaceholders(item, caseDir, configDir));
        return core::Value(std::move(out));
    }
    if (const auto* t = v.get<core::Table>()) {
        core::Table out;
        for (const auto& [k, item] : *t) out.emplace_back(k, withPlaceholders(item, caseDir, configDir));
        return core::Value(std::move(out));
    }
    return v;
}

std::string compareList(const std::string& what, const core::Value* expected, const core::Array& actual,
                        const fs::path& caseDir, const fs::path& configDir) {
    const core::Array none;
    const core::Value want = masked(core::Value(expected && expected->isArray() ? *expected->get<core::Array>() : none));
    const core::Value got = withPlaceholders(core::Value(actual), caseDir, configDir);
    if (sameValue(want, got)) return {};
    return what + ": expected " + core::toJson(want, -1) + ", got " + core::toJson(got, -1) + "\n";
}

}  // namespace

bool sameValue(const core::Value& a, const core::Value& b) {
    if (isNumber(a) && isNumber(b)) {
        const double x = numberOf(a);
        const double y = numberOf(b);
        if (std::isnan(x) || std::isnan(y)) return std::isnan(x) && std::isnan(y);
        if (x == y) return true;
        return std::abs(x - y) <= 1e-9 * std::max(std::abs(x), std::abs(y));
    }
    if (a.type() != b.type()) return false;
    if (const auto* x = a.get<core::Array>()) {
        const auto& y = *b.get<core::Array>();
        if (x->size() != y.size()) return false;
        for (std::size_t k = 0; k < x->size(); ++k) {
            if (!sameValue((*x)[k], y[k])) return false;
        }
        return true;
    }
    if (const auto* x = a.get<core::Table>()) {
        const auto& y = *b.get<core::Table>();
        if (x->size() != y.size()) return false;
        for (const auto& [key, value] : *x) {
            const core::Value* other = b.find(key);
            if (!other || !sameValue(value, *other)) return false;
        }
        return true;
    }
    return a == b;
}

std::string checkOutcome(const GuiCase& c, const CaseRunner& runner, const std::vector<std::string>& log,
                         const fs::path& dataDir, const fs::path& caseDir, const fs::path& configDir) {
    std::string problems;
    const core::Value& outcome = c.outcome;
    if (!outcome.isTable()) {  // a case of the harness's own format: its files only
        return problems;
    }
    const auto entry = [&](const char* key) { return outcome.isTable() ? outcome.find(key) : nullptr; };
    if (c.log) {
        std::vector<std::string> got;
        for (const std::string& line : log) got.push_back(placeholders(line, caseDir, configDir));
        std::vector<std::string> want;
        for (const std::string& line : *c.log) want.push_back(maskParserWords(line));
        if (got != want) {
            std::string text = "log: expected\n";
            for (const std::string& line : want) text += "    " + line + "\n";
            text += "  got\n";
            for (const std::string& line : got) text += "    " + line + "\n";
            problems += text;
        }
    }
    problems += compareList("questions", entry("questions"), runner.prompter().questions, caseDir, configDir);
    problems += compareList("dialogs", entry("dialogs"), runner.prompter().dialogs, caseDir, configDir);
    problems += compareList("snapshots", entry("snapshots"), runner.snapshots(), caseDir, configDir);
    {
        core::Array lines;
        for (const std::string& line : runner.prompter().stderrLines) lines.push_back(core::Value(line));
        problems += compareList("stderr", entry("stderr"), lines, caseDir, configDir);
    }
    const fs::path settings = core::userSettingsPath();
    const core::Value* expected = entry("settings");
    std::error_code ec;
    if (expected && expected->isString()) {
        if (!fs::exists(settings, ec)) {
            problems += "settings: " + core::pathToUtf8(settings) + " was not written\n";
        } else {
            // the comments line for line and the values; arrays wrap by the length of the paths in them,
            // which differ between the generator's folders and these; the placeholders go into the
            // parsed values (in the text, a TOML string escapes a Windows path's '\\')
            const std::string want = normalise(readFile(dataDir / core::pathFromUtf8(expected->stringOr())), {});
            const std::string written = readFile(settings);
            const std::string got = placeholders(written, caseDir, configDir);
            const auto comments = [](const std::string& text) {
                std::istringstream in(text);
                std::string out;
                for (std::string line; std::getline(in, line);) {
                    if (!line.empty() && line[0] == '#') out += line + "\n";
                }
                return out;
            };
            try {
                if (const std::string diff = compareText(comments(got), comments(want)); !diff.empty()) {
                    problems += "settings (comments): " + diff + "\n";
                } else if (!sameValue(withPlaceholders(core::parseToml(written, "settings.toml"), caseDir, configDir),
                                      core::parseToml(want, "expected"))) {
                    problems += "settings: expected\n" + want + "got\n" + got;
                }
            } catch (const std::exception& e) {
                problems += std::string("settings: ") + e.what() + "\n";
            }
        }
    } else if (fs::exists(settings, ec)) {
        problems += "settings: " + core::pathToUtf8(settings) + " was written, but the case expects none\n";
    }
    return problems;
}

}  // namespace mdbridge::guitest
