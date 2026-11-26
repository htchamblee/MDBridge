#include "gui_harness.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>

#include <QAbstractButton>
#include <QAbstractSpinBox>
#include <QApplication>
#include <QComboBox>
#include <QImage>
#include <QLineEdit>
#include <QPixmap>
#include <QPlainTextEdit>
#include <QTableView>
#include <QTest>
#include <QTextEdit>
#include <QWidget>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/main_window.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/project_document.hpp"
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

// The text a user types for `value`.
QString typedText(const core::Value& value, const char* widgetType) {
    if (const auto* s = value.get<std::string>()) return toQString(*s);
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
    c.actions = *actions->get<core::Array>();
    if (const core::Value* start = root.find("start")) {
        if (!start->isTable() || !start->find("project") || !start->find("project")->isString()) {
            throw CaseError(where + ": 'start' must be {\"project\": name}");
        }
        c.startProject = start->find("project")->stringOr();
    }
    if (const core::Value* outcome = root.find("outcome")) {
        c.outcome = *outcome;
    }
    for (std::size_t k = 0; k < c.actions.size(); ++k) {
        if (!c.actions[k].isTable() || !c.actions[k].find("action") || !c.actions[k].find("action")->isString()) {
            throw CaseError(where + ": action " + std::to_string(k + 1) + " must be an object with an \"action\" string");
        }
    }
    return c;
}

std::vector<std::string> unsupportedNeeds(const GuiCase& c, const std::vector<std::string>& extra) {
    // grows as the window and this runner do: "type" and "table" are run (their cases still need the
    // settings pages of wave 1, and "start" / "lossless" / "snapshot", so they stay skipped)
    static const std::vector<std::string> supported{"table", "type"};
    std::vector<std::string> out;
    for (const std::string& need : c.needs) {
        if (std::find(supported.begin(), supported.end(), need) == supported.end() &&
            std::find(extra.begin(), extra.end(), need) == extra.end()) {
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

fs::path screenshotDir(const fs::path& fallback) {
    if (const char* dir = std::getenv("MDBRIDGE_GUI_SCREENSHOTS"); dir && *dir) {
        return core::pathFromUtf8(dir);
    }
    return fallback;
}

std::string normalise(std::string text, const fs::path& caseDir) {
    replaceAll(text, "\r\n", "\n");
    if (!caseDir.empty()) {
        const fs::path dir = caseDir.lexically_normal();
        replaceAll(text, genericUtf8(dir), "{case_dir}");
        replaceAll(text, core::pathToUtf8(dir), "{case_dir}");
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

namespace {

// The main window as a case host.
class WindowHost : public CaseHost {
public:
    explicit WindowHost(app::MainWindow& window) : m_window(window) {}
    QWidget* widget() override { return &m_window; }
    bool openProject(const fs::path& path) override { return m_window.openProject(path, false); }
    bool saveProject(const fs::path& path) override { return m_window.saveProject(path, false); }
    bool settled() const override { return !m_window.loading() && m_window.preview()->settled(); }
    core::ProjectDocument* document() override { return m_window.document(); }

private:
    app::MainWindow& m_window;
};

// The value as compared: tables by key (order free), arrays item by item, numbers by type.
bool sameValue(const core::Value& a, const core::Value& b) {
    if (a.type() != b.type()) {
        return false;
    }
    if (const core::Table* ta = a.get<core::Table>()) {
        const core::Table& tb = *b.get<core::Table>();
        if (ta->size() != tb.size()) return false;
        for (const auto& [key, value] : *ta) {
            const core::Value* other = b.find(key);
            if (!other || !sameValue(value, *other)) return false;
        }
        return true;
    }
    if (const core::Array* xa = a.get<core::Array>()) {
        const core::Array& xb = *b.get<core::Array>();
        if (xa->size() != xb.size()) return false;
        for (std::size_t k = 0; k < xa->size(); ++k) {
            if (!sameValue((*xa)[k], xb[k])) return false;
        }
        return true;
    }
    if (const double* da = a.get<double>()) {
        const double db = *b.get<double>();
        return *da == db || (std::isnan(*da) && std::isnan(db));
    }
    return a == b;
}

}  // namespace

CaseRunner::CaseRunner(app::MainWindow& window, fs::path caseDir, fs::path screenshots, std::string caseName)
    : CaseRunner(*new WindowHost(window), std::move(caseDir), std::move(screenshots), std::move(caseName)) {
    m_ownHost.reset(&m_host);
}

CaseRunner::CaseRunner(CaseHost& host, fs::path caseDir, fs::path screenshots, std::string caseName)
    : m_host(host), m_caseDir(std::move(caseDir)), m_screenshotDir(std::move(screenshots)),
      m_caseName(std::move(caseName)) {
    std::replace(m_caseName.begin(), m_caseName.end(), '/', '_');
}

CaseRunner::~CaseRunner() = default;

void CaseRunner::start(const GuiCase& c) {
    if (c.startProject.empty()) {
        return;
    }
    const fs::path file = m_caseDir / core::pathFromUtf8(c.startProject);
    if (!m_host.openProject(file)) {
        throw CaseError("start: cannot open " + core::pathToUtf8(file));
    }
}

void CaseRunner::snapshot(const core::Value& action) {
    m_host.commitPendingEdits();
    core::Table record{{"action", core::Value(static_cast<std::int64_t>(m_action))}};
    if (const core::Value* keys = action.find("keys"); keys && keys->isArray() && !keys->get<core::Array>()->empty()) {
        core::ProjectDocument* doc = m_host.document();
        if (!doc) {
            throw CaseError("no project is open");
        }
        core::Table values;
        for (const core::Value& key : *keys->get<core::Array>()) {
            values.emplace_back(key.stringOr(), doc->get(key.stringOr()));
        }
        record.emplace_back("values", core::Value(std::move(values)));
    }
    if (const core::Value* page = action.find("page"); page && page->get<bool>() && *page->get<bool>()) {
        record.emplace_back("page", core::Value(m_host.currentPage().toStdString()));
    }
    if (const core::Value* modified = action.find("modified"); modified && modified->get<bool>() && *modified->get<bool>()) {
        record.emplace_back("modified", core::Value(m_host.document() && m_host.document()->isModified()));
    }
    m_snapshots.push_back({m_action, core::Value(std::move(record))});
}

std::string compareSnapshots(const GuiCase& c, const std::vector<Snapshot>& recorded, std::vector<std::string>* skipped) {
    const core::Value* expected = c.outcome.isTable() ? c.outcome.find("snapshots") : nullptr;
    const core::Array none;
    const core::Array& want = expected && expected->isArray() ? *expected->get<core::Array>() : none;
    if (want.size() != recorded.size()) {
        return "outcome.snapshots: expected " + std::to_string(want.size()) + " snapshot(s), recorded " +
               std::to_string(recorded.size()) + "\n";
    }
    std::string problems;
    for (std::size_t k = 0; k < want.size(); ++k) {
        const core::Value& w = want[k];
        const core::Value& got = recorded[k].record;
        const std::string where = "snapshot of action " + std::to_string(recorded[k].action);
        if (const core::Value* table = w.isTable() ? &w : nullptr) {
            for (const auto& [key, value] : *table->get<core::Table>()) {
                if (key == "action") {
                    continue;
                }
                if (key == "slots" || key == "density_status") {
                    if (skipped) skipped->push_back(where + ": " + key + " (no layout / structure page in this host)");
                    continue;
                }
                const core::Value* mine = got.find(key);
                if (key == "values" && mine) {
                    for (const auto& [field, v] : *value.get<core::Table>()) {
                        const core::Value* m = mine->find(field);
                        if (!m || !sameValue(*m, v)) {
                            problems += where + ": " + field + ": expected " + core::toJson(v, -1) + ", got " +
                                        (m ? core::toJson(*m, -1) : std::string("nothing")) + "\n";
                        }
                    }
                } else if (!mine || !sameValue(*mine, value)) {
                    problems += where + ": " + key + ": expected " + core::toJson(value, -1) + ", got " +
                                (mine ? core::toJson(*mine, -1) : std::string("nothing")) + "\n";
                }
            }
        }
    }
    return problems;
}

void CaseRunner::run(const core::Array& actions) {
    for (std::size_t k = 0; k < actions.size(); ++k) {
        try {
            m_action = k;
            runAction(actions[k]);
        } catch (const std::exception& e) {
            throw CaseError("action " + std::to_string(k + 1) + " (" + actions[k].find("action")->stringOr() +
                            "): " + e.what());
        }
    }
}

void CaseRunner::waitUntilLoaded(int timeoutMs) {
    const bool done = QTest::qWaitFor([this] { return m_host.settled(); }, timeoutMs);
    if (!done) {
        throw CaseError("the project was not loaded in " + std::to_string(timeoutMs / 1000) + " s");
    }
}

void CaseRunner::runAction(const core::Value& action) {
    const std::string& what = stringEntry(action, "action", "action");
    if (what == "open") {
        const fs::path file = m_caseDir / core::pathFromUtf8(stringEntry(action, "path", what));
        if (!m_host.openProject(file)) {
            throw CaseError("cannot open " + core::pathToUtf8(file));
        }
    } else if (what == "wait") {
        int timeoutMs = kDefaultWaitMs;
        if (const core::Value* t = action.find("timeout_s")) {
            if (const auto* i = t->get<std::int64_t>()) timeoutMs = static_cast<int>(*i * 1000);
            else if (const auto* d = t->get<double>()) timeoutMs = static_cast<int>(*d * 1000.0);
            else throw CaseError("'timeout_s' must be a number");
        }
        waitUntilLoaded(timeoutMs);
    } else if (what == "set") {
        const core::Value* value = action.find("value");
        if (!value) {
            throw CaseError("'value' is missing");
        }
        m_host.reveal(toQString(stringEntry(action, "key", what)));
        setKey(m_host.widget(), toQString(stringEntry(action, "key", what)), *value);
    } else if (what == "type") {
        m_host.reveal(toQString(stringEntry(action, "key", what)));
        typeKey(m_host.widget(), toQString(stringEntry(action, "key", what)), toQString(stringEntry(action, "text", what)));
    } else if (what == "table") {
        const auto integer = [&](std::string_view k, int fallback) {
            const core::Value* v = action.find(k);
            if (!v) return fallback;
            const auto* i = v->get<std::int64_t>();
            if (!i) throw CaseError("'" + std::string(k) + "' must be a whole number");
            return static_cast<int>(*i);
        };
        m_host.reveal(toQString(stringEntry(action, "key", what)));
        tableAction(m_host.widget(), toQString(stringEntry(action, "key", what)), stringEntry(action, "do", what),
                    integer("row", -1), integer("by", 0));
    } else if (what == "click") {
        m_host.reveal(toQString(stringEntry(action, "key", what)));
        QWidget* widget = findField(m_host.widget(), toQString(stringEntry(action, "key", what)));
        auto* button = qobject_cast<QAbstractButton*>(widget);
        if (!button) {
            throw CaseError(widget->objectName().toStdString() + " is no button");
        }
        QTest::mouseClick(button, Qt::LeftButton);
    } else if (what == "save") {
        const std::string path = optionalString(action, "path");
        const fs::path file = path.empty() ? fs::path() : m_caseDir / core::pathFromUtf8(path);
        m_host.commitPendingEdits();
        if (!m_host.saveProject(file)) {
            throw CaseError("the save failed (see the output above)");
        }
    } else if (what == "screenshot") {
        const std::string widgetName = optionalString(action, "widget");
        QWidget* widget = widgetName.empty() ? m_host.widget() : findField(m_host.widget(), toQString(widgetName));
        m_shots.push_back(saveScreenshot(widget, m_screenshotDir, m_caseName + "-" + stringEntry(action, "name", what)));
    } else if (what == "snapshot") {
        snapshot(action);
    } else {
        throw CaseError("unknown action (tests/data/README.md, \"gui\", lists the actions)");
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

}  // namespace mdbridge::guitest
