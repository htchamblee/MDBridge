// The GUI's bound field widgets (src/app/gui/forms, milestone 5 W0-C), Qt Test on the offscreen
// platform: a form of every field of every fixture project saves the file back byte for byte
// without an edit; edits store one undo step and undo / redo refresh the widgets; AP-01 .. AP-05,
// AP-10, AP-15 (docs/specs/m5-plan.md section 4); the mouse wheel (MDATK 0.0.4
// test_value_widgets_ignore_the_mouse_wheel); the inputs the GUI harness uses (typing, picking
// items by project value, Space on check boxes, table cells).

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <QAbstractSlider>
#include <QApplication>
#include <QDir>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QTableView>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <QtTest>

#include "gui/forms/field_factory.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/forms/list_editor.hpp"
#include "mdbridge/core/clock.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project_io.hpp"
#include "mdbridge/core/project_file.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::app::forms;
using mdbridge::test::readFile;
using mdbridge::test::TempDir;
namespace fs = std::filesystem;

namespace {

const fs::path kData = MDBRIDGE_TEST_DATA;
const fs::path kProjects = kData / "projects";

std::string unixLines(const std::string& text) {
    std::string out;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (!(text[i] == '\r' && i + 1 < text.size() && text[i + 1] == '\n')) {
            out += text[i];
        }
    }
    return out;
}

class FixedClock {
public:
    explicit FixedClock(core::DateTime now) : m_previous(core::setClock(core::Clock{[now] { return now; }, {}})) {}
    ~FixedClock() { core::setClock(m_previous); }
    FixedClock(const FixedClock&) = delete;
    FixedClock& operator=(const FixedClock&) = delete;

private:
    core::Clock m_previous;
};

core::DateTime savedStamp(const fs::path& path) {
    const core::Value* saved = core::readProject(path).header.find("saved");
    if (saved && saved->get<core::DateTime>()) {
        return *saved->get<core::DateTime>();
    }
    return mdbridge::test::localDateTime(2026, 11, 25, 12, 0, 0);  // a file without a stamp
}

// The project files MDATK wrote (folder, file), as tests/core/test_project_document.cpp lists them.
std::vector<std::pair<fs::path, std::string>> fixtureProjects() {
    std::vector<std::pair<fs::path, std::string>> out;
    for (const char* name : {"default.mdatk.toml", "full.mdatk.toml", "split.mdatk.toml", "handwritten.mdatk.toml",
                             "datetimes.mdatk.toml"}) {
        out.emplace_back(kProjects, name);
    }
    const fs::path golden = kData / "project-file" / "golden";
    for (const char* name : {"default.mdatk.toml", "rich.mdatk.toml", "rich.resaved.mdatk.toml",
                             "rich.with-log.mdatk.toml", "split-nolog.after.mdatk.toml", "split.after.mdatk.toml",
                             "split.before.mdatk.toml", "split.resaved.mdatk.toml", "truncation.mdatk.toml"}) {
        out.emplace_back(golden, name);
    }
    const fs::path workflow = kData / "project-file" / "workflow";
    for (const auto& entry : fs::directory_iterator(workflow)) {
        const std::string name = entry.path().filename().string();
        if (name.size() > 11 && name.substr(name.size() - 11) == ".mdatk.toml") {
            out.emplace_back(workflow, name);
        }
    }
    for (const char* name : {"plan.mdatk.toml", "sidecar.mdatk.toml"}) {
        out.emplace_back(kData / "project-file" / "mdatk-0.1.1", name);
    }
    return out;
}

void copyFolder(const fs::path& from, const fs::path& to) {
    for (const auto& entry : fs::directory_iterator(from)) {
        if (entry.is_regular_file()) {
            fs::copy_file(entry.path(), to / entry.path().filename(), fs::copy_options::overwrite_existing);
        }
    }
}

// A window with a widget for every field of the document (formFieldPaths), shown and active.
struct Form {
    explicit Form(core::ProjectDocument& doc) : context(doc) {
        window = std::make_unique<QWidget>();
        auto* outer = new QVBoxLayout(window.get());
        auto* scroll = new QScrollArea(window.get());
        scroll->setWidgetResizable(true);
        auto* page = new QWidget;
        auto* layout = new QVBoxLayout(page);
        for (const std::string& path : formFieldPaths(doc.project())) {
            if (FieldWidget* field = createFieldWidget(context, path, page)) {
                layout->addWidget(field);
                fields.push_back(field);
            }
        }
        scroll->setWidget(page);
        outer->addWidget(scroll);
        window->resize(600, 800);
        window->show();
        window->activateWindow();
        active = QTest::qWaitForWindowActive(window.get());
        context.refreshProblems();
    }

    template <class T = FieldWidget>
    T* field(const std::string& path) const {
        for (FieldWidget* f : fields) {
            if (f->path() == path) {
                return qobject_cast<T*>(f);
            }
        }
        return nullptr;
    }

    FormContext context;
    std::unique_ptr<QWidget> window;
    std::vector<FieldWidget*> fields;
    bool active = false;
};

// One field alone in a shown window.
struct One {
    One(core::ProjectDocument& doc, const std::string& path) : context(doc) {
        window = std::make_unique<QWidget>();
        auto* layout = new QVBoxLayout(window.get());
        field = createFieldWidget(context, path, window.get());
        other = new QLineEdit(window.get());  // somewhere else to click
        if (field) {
            layout->addWidget(field);
        }
        layout->addWidget(other);
        window->show();
        window->activateWindow();
        active = QTest::qWaitForWindowActive(window.get());
        context.refreshProblems();
    }
    FormContext context;
    std::unique_ptr<QWidget> window;
    FieldWidget* field = nullptr;
    QLineEdit* other = nullptr;
    bool active = false;
};

std::size_t undoSteps(core::ProjectDocument& doc) {
    std::size_t n = 0;
    while (doc.canUndo()) {
        doc.undo();
        ++n;
    }
    for (std::size_t k = 0; k < n; ++k) {
        doc.redo();
    }
    return n;
}

// As the GUI harness types into a line edit, spin box or editable combo: select all, type, Enter.
void typeLine(QWidget* widget, const QString& text) {
    if (auto* spin = qobject_cast<QAbstractSpinBox*>(widget)) {
        spin->selectAll();
    } else if (auto* combo = qobject_cast<QComboBox*>(widget)) {
        combo->lineEdit()->selectAll();
        widget = combo->lineEdit();
    } else if (auto* edit = qobject_cast<QLineEdit*>(widget)) {
        edit->selectAll();
    }
    if (text.isEmpty()) {
        QTest::keyClick(widget, Qt::Key_Backspace);
    }
    QTest::keyClicks(widget, text);
    QTest::keyClick(widget, Qt::Key_Return);
}

std::string text(const core::ProjectDocument& doc, const std::string& path) { return doc.text(path); }

}  // namespace

class FormsTest : public QObject {
    Q_OBJECT

private slots:
    void everyFieldOfEveryFixtureSavesBackByteForByte();
    void numbersAreNeverClamped();
    void numberArrowsStepWithinTheRangeHint();
    void optionalNumbersKeepEveryDigit();
    void coloursAreKeptAsTyped();
    void chartFormatsAreKeptAsAList();
    void unknownChoicesAreKeptAndShown();
    void formErrorsAreLocatedAndKeepTheValue();
    void documentProblemsFlagTheirField();
    void oneEditIsOneUndoStepAndUndoRefreshes();
    void multiLineTextCommitsOnLeaving();
    void checkBoxesAndEditableCombos();
    void pathFieldsBrowseWithMdatksFilters();
    void listEditorRowsAndUndo();
    void tableCellsCommitOnFocusLossAndEnter();
    void valueWidgetsIgnoreTheMouseWheel();
    void everyEditorIsNamedAndReachable();
};

// Load -> save is the identity for every field: a form of every field of each fixture project,
// each field given Enter (and the pending edits committed) without typing, saves the file back
// byte for byte, with no undo step; each widget shows the document's own text.
void FormsTest::everyFieldOfEveryFixtureSavesBackByteForByte() {
    std::size_t files = 0;
    for (const auto& [folder, name] : fixtureProjects()) {
        TempDir dir;
        copyFolder(folder, dir.path());
        const fs::path path = dir / name;
        const std::string before = readFile(path);
        core::ProjectDocument doc;
        doc.open(path);
        Form form(doc);
        QVERIFY2(form.fields.size() > 150, name.c_str());
        for (FieldWidget* field : form.fields) {
            if (qobject_cast<ListEditor*>(field)) {
                QCOMPARE(qobject_cast<ListEditor*>(field)->rowCount(),
                         static_cast<int>(doc.get(field->path()).size()));
                continue;
            }
            QString shown = toQString(doc.text(field->path()));
            if (qobject_cast<MultiLineField*>(field)) {
                // a text edit shows a lone carriage return as a line break (kept unless edited)
                shown.replace(QStringLiteral("\r\n"), QStringLiteral("\n")).replace(QLatin1Char('\r'), QLatin1Char('\n'));
            }
            QCOMPARE(field->displayedText(), shown);
            if (qobject_cast<QLineEdit*>(field->editor()) || qobject_cast<QAbstractSpinBox*>(field->editor())) {
                field->editor()->setFocus();
                QTest::keyClick(field->editor(), Qt::Key_Return);
            }
        }
        form.context.commitPendingEdits();
        QVERIFY2(!doc.canUndo(), name.c_str());
        QVERIFY2(!doc.isModified(), name.c_str());
        fs::copy_file(path, dir / "original.toml");
        {
            const FixedClock clock(savedStamp(path));
            doc.save();
        }
        std::string expected = unixLines(before);
        if (name == "split.after.mdatk.toml" || name == "sidecar.mdatk.toml") {
            expected.pop_back();  // MDATK's log move leaves a blank line that its next save drops too
        }
        if (name == "handwritten.mdatk.toml" || name == "datetimes.mdatk.toml") {
            // not written by MDATK: the save is MDATK's form of the same values (W0-B's test)
            core::ProjectDocument again;
            again.open(path);
            QCOMPARE(toQString(core::toJson(core::toPlain(again.resolvedProject()), 1)),
                     toQString(core::toJson(core::toPlain(core::loadProject(dir / "original.toml")), 1)));
        } else {
            QCOMPARE(toQString(unixLines(readFile(path))), toQString(expected));
        }
        ++files;
    }
    QVERIFY(files >= 20);
}

// AP-01: values outside MDATK's spin ranges are stored as typed and flagged, never clamped.
void FormsTest::numbersAreNeverClamped() {
    core::ProjectDocument doc;
    Form form(doc);
    for (const auto& [path, typed, stored] :
         std::vector<std::tuple<std::string, QString, std::string>>{{"analyses.contacts.cutoff", "25", "25.0"},
                                                                    {"analyses.hbonds.angle_cutoff", "60", "60.0"},
                                                                    {"output.video.width", "7682", "7682"},
                                                                    {"visual.zoom", "0.01", "0.01"},
                                                                    {"analyses.energy.cutoff", "60", "60.0"}}) {
        auto* field = form.field<NumberField>(path);
        QVERIFY(field);
        typeLine(field->editor(), typed);
        QCOMPARE(text(doc, path), stored);
        QCOMPARE(field->displayedText(), toQString(stored));
        if (form.context.problemsAt(path).empty()) {  // validate() allows it: the range hint flags it
            QCOMPARE(field->flag().second, QStringLiteral("warning"));
            QVERIFY(field->flag().first.contains(QStringLiteral("outside the usual range")));
        } else {  // validate()'s own limit
            QCOMPARE(field->flag().first, toQString(form.context.problemsAt(path).front().message));
        }
        QVERIFY(field->editor()->toolTip().contains(QStringLiteral("Usual range")));
    }
    QCOMPARE(undoSteps(doc), std::size_t{5});
    // inside the range: no flag
    auto* cutoff = form.field<NumberField>("analyses.contacts.cutoff");
    typeLine(cutoff->editor(), "4.5");
    QCOMPARE(cutoff->flag().first, QString());
    // a value read from a file is shown as read and stays unless edited
    doc.set("analyses.contacts.cutoff", core::Value(25.0));
    QCOMPARE(cutoff->displayedText(), QStringLiteral("25.0"));
}

void FormsTest::numberArrowsStepWithinTheRangeHint() {
    core::ProjectDocument doc;
    One one(doc, "analyses.contacts.cutoff");
    auto* field = qobject_cast<NumberField*>(one.field);
    QVERIFY(field);
    doc.set("analyses.contacts.cutoff", core::Value(4.5));
    field->box()->stepBy(1);
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("4.6"));  // 2 decimals: no 4.6000000000000005
    field->box()->stepBy(-3);
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("4.3"));
    doc.set("analyses.contacts.cutoff", core::Value(19.95));
    field->box()->stepBy(1);
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("20.0"));  // stops at the hint
    QVERIFY(!(field->box()->stepEnabled() & QAbstractSpinBox::StepUpEnabled));
    doc.set("analyses.contacts.cutoff", core::Value(25.0));  // outside: kept; down steps from it
    QVERIFY(!(field->box()->stepEnabled() & QAbstractSpinBox::StepUpEnabled));
    field->box()->stepBy(-1);
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("24.9"));
    // integers
    One width(doc, "output.video.width");
    auto* w = qobject_cast<NumberField*>(width.field);
    doc.set("output.video.width", core::Value(std::int64_t{1279}));
    w->box()->stepBy(1);
    QCOMPARE(text(doc, "output.video.width"), std::string("1281"));
    // text the field cannot hold: the arrows start from the hint's minimum, no crash
    typeLine(w->editor(), "abc");
    QCOMPARE(text(doc, "output.video.width"), std::string("1281"));
}

// AP-02: optional numbers are shown and stored with every digit (no %g).
void FormsTest::optionalNumbersKeepEveryDigit() {
    core::ProjectDocument doc;
    Form form(doc);
    typeLine(form.field("inputs.dt_ps")->editor(), "0.0123456789");
    typeLine(form.field("analyses.energy.target_value")->editor(), "-1234567.5");
    QCOMPARE(text(doc, "inputs.dt_ps"), std::string("0.0123456789"));
    QCOMPARE(text(doc, "analyses.energy.target_value"), std::string("-1234567.5"));
    QCOMPARE(form.field("inputs.dt_ps")->displayedText(), QStringLiteral("0.0123456789"));
    // "" unsets; the placeholder is shown
    typeLine(form.field("inputs.dt_ps")->editor(), "");
    QVERIFY(doc.get("inputs.dt_ps").isNull());
    QCOMPARE(qobject_cast<QLineEdit*>(form.field("inputs.stop")->editor())->placeholderText(), QStringLiteral("last"));
}

// AP-03: any Matplotlib colour text is kept as typed; the swatch shows it; the picker writes
// #rrggbb; text Matplotlib cannot read is flagged with its message.
void FormsTest::coloursAreKeptAsTyped() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* field = form.field<ColorField>("visual.base_color");
    QVERIFY(field);
    for (const QString& c : {QStringLiteral("tab:orange"), QStringLiteral("#ff000080"), QStringLiteral("gold"),
                             QStringLiteral("0.5"), QStringLiteral("C3")}) {
        typeLine(field->editor(), c);
        QCOMPARE(toQString(text(doc, "visual.base_color")), c);
        QVERIFY(field->swatchColor().isValid());
        QCOMPARE(field->flag().first, QString());
    }
    QCOMPARE(field->swatchColor().alpha(), 255);
    typeLine(field->editor(), "#ff000080");
    QCOMPARE(field->swatchColor().alpha(), 128);
    typeLine(field->editor(), "notacolour");
    QCOMPARE(text(doc, "visual.base_color"), std::string("notacolour"));
    QVERIFY(!field->swatchColor().isValid());
    QCOMPARE(field->flag(), std::make_pair(QStringLiteral("Invalid RGBA argument: 'notacolour'"), QStringLiteral("error")));
    // the picker
    QColor initial;
    const auto previous = setColorDialogFunction([&](QWidget*, const QColor& start, const QString&) {
        initial = start;
        return QColor(0x12, 0x34, 0x56);
    });
    typeLine(field->editor(), "tab:orange");
    field->pickButton()->click();
    QCOMPARE(initial, matplotlibColor(QStringLiteral("tab:orange")));
    QCOMPARE(text(doc, "visual.base_color"), std::string("#123456"));
    setColorDialogFunction([](QWidget*, const QColor&, const QString&) { return QColor(); });  // cancelled
    field->pickButton()->click();
    QCOMPARE(text(doc, "visual.base_color"), std::string("#123456"));
    setColorDialogFunction(previous);
    // a colour that may be empty
    auto* charts = form.field<ColorField>("output.charts.color");
    QVERIFY(charts);
    typeLine(charts->editor(), "");
    QCOMPARE(charts->flag().first, QString());
    QCOMPARE(qobject_cast<QLineEdit*>(charts->editor())->placeholderText(), QStringLiteral("palette"));
}

// AP-04: the chart formats are a list kept as it is (a page may show check boxes; the list widget
// keeps entries it has no box for).
void FormsTest::chartFormatsAreKeptAsAList() {
    core::ProjectDocument doc;
    Form form(doc);
    typeLine(form.field("output.charts.formats")->editor(), "jpg, svg");
    QCOMPARE(text(doc, "output.charts.formats"), std::string("jpg, svg"));
    typeLine(form.field("output.charts.formats")->editor(), "");
    QCOMPARE(doc.get("output.charts.formats").size(), std::size_t{0});
}

// AP-05: a value not among the choices stays, shown as "<value> (unknown)" and flagged; never
// replaced (also not by the choice of a project opened before).
void FormsTest::unknownChoicesAreKeptAndShown() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* field = form.field<ComboField>("visual.renderer");
    QVERIFY(field);
    QComboBox* combo = field->comboBox();
    // items carry the project value as data (the harness picks by it)
    QVERIFY(combo->findData(QStringLiteral("software")) >= 0);
    QCOMPARE(combo->currentData().toString(), QStringLiteral("auto"));
    const int software = combo->findData(QStringLiteral("software"));
    combo->setCurrentIndex(software);
    emit combo->activated(software);
    QCOMPARE(text(doc, "visual.renderer"), std::string("software"));
    doc.set("visual.renderer", core::Value(std::string("opengl")));
    form.context.refreshProblems();
    QCOMPARE(combo->currentText(), QStringLiteral("opengl (unknown)"));
    QCOMPARE(field->displayedText(), QStringLiteral("opengl"));
    QCOMPARE(field->flag().second, QStringLiteral("error"));  // validate() reports it (RN-20)
    QVERIFY(field->flag().first.contains(QStringLiteral("opengl")));
    // a combo of a field validate() does not check flags it itself, as a warning
    auto* engine = form.field<ComboField>("analyses.basepairs.engine");
    doc.set("analyses.basepairs.engine", core::Value(std::string("barnaba")));
    form.context.refreshProblems();
    QCOMPARE(engine->comboBox()->currentText(), QStringLiteral("barnaba (unknown)"));
    QVERIFY(!engine->flag().first.isEmpty());
    // undo / redo show the values again
    doc.undo();
    QCOMPARE(engine->comboBox()->currentData().toString(), QStringLiteral("builtin"));
    QCOMPARE(engine->comboBox()->count(), 1);
    doc.redo();
    QCOMPARE(engine->displayedText(), QStringLiteral("barnaba"));
    // a project opened after another keeps its own value
    TempDir dir;
    core::ProjectDocument other;
    other.set("visual.renderer", core::Value(std::string("vulkan2")));
    other.saveAs(dir / "a.mdatk.toml");
    doc.open(dir / "a.mdatk.toml");
    QCOMPARE(field->displayedText(), QStringLiteral("vulkan2"));
    form.context.commitPendingEdits();
    QVERIFY(!doc.isModified());
}

// AP-10: text a field cannot hold is refused with the located message, shown with the field; the
// document keeps its value; Esc shows the value again.
void FormsTest::formErrorsAreLocatedAndKeepTheValue() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* stop = form.field<LineField>("inputs.stop");
    typeLine(stop->editor(), "x");
    QVERIFY(doc.get("inputs.stop").isNull());
    QVERIFY(!doc.canUndo());
    QCOMPARE(stop->displayedText(), QStringLiteral("x"));
    QCOMPARE(stop->flag(), std::make_pair(QStringLiteral("'x' is not a whole number"), QStringLiteral("error")));
    QCOMPARE(form.context.formErrors().at("inputs.stop").text(),
             std::string("Input & cleaning › Frames › Stop: 'x' is not a whole number"));
    QVERIFY(stop->editor()->toolTip().startsWith(QStringLiteral("'x' is not a whole number")));
    // Esc: the stored value again, no error
    stop->editor()->setFocus();
    QTest::keyClicks(stop->editor(), "y");
    QTest::keyClick(stop->editor(), Qt::Key_Escape);
    QCOMPARE(stop->displayedText(), QString());
    QVERIFY(form.context.formErrors().empty());
    QCOMPARE(stop->flag().first, QString());
    // a number field and a JSON field
    typeLine(form.field("analyses.contacts.cutoff")->editor(), "4,5");
    QCOMPARE(form.field("analyses.contacts.cutoff")->flag().first, QStringLiteral("'4,5' is not a number"));
    typeLine(form.field("analyses.contacts.cutoff")->editor(), "4.5");
    QCOMPARE(form.field("analyses.contacts.cutoff")->flag().first, QString());
    typeLine(form.field("layout.panels.0.options")->editor(), "[1, 2]");
    QCOMPARE(form.field("layout.panels.0.options")->flag().first, QStringLiteral("options must be a JSON object"));
    QCOMPARE(text(doc, "layout.panels.0.options"), std::string());
    QCOMPARE(undoSteps(doc), std::size_t{0});  // refused text is no step
}

// Problems of the document flag their field (after the problems are recomputed: once per change).
void FormsTest::documentProblemsFlagTheirField() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* group = form.field<EditableComboField>("analyses.contacts.group_a");
    QVERIFY(group);
    typeLine(group->editor(), "Nope");
    QCOMPARE(text(doc, "analyses.contacts.group_a"), std::string("Nope"));
    QTRY_VERIFY(!group->flag().first.isEmpty());
    QCOMPARE(group->flag().second, QStringLiteral("error"));
    QVERIFY(group->flag().first.contains(QStringLiteral("Nope")));
    doc.undo();
    QTRY_COMPARE(group->flag().first, QString());
    // a problem of a record field tints its cell in the list editor
    auto* groups = form.field<ListEditor>("groups");
    QVERIFY(groups);
    doc.setText("groups.1.name", "Aptamer");
    form.context.refreshProblems();
    const QModelIndex cell = groups->cellIndex("groups.1.name");
    QVERIFY(cell.isValid());
    QVERIFY(!cell.data(Qt::ToolTipRole).toString().isEmpty());
    QVERIFY(cell.data(Qt::BackgroundRole).isValid());
}

// One user edit is one undo step (text stored when editing finishes, not per key); undo and redo
// refresh the widgets without creating steps.
void FormsTest::oneEditIsOneUndoStepAndUndoRefreshes() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* cutoff = form.field("analyses.contacts.cutoff");
    cutoff->editor()->setFocus();
    qobject_cast<QAbstractSpinBox*>(cutoff->editor())->selectAll();
    QTest::keyClicks(cutoff->editor(), "4.75");
    QVERIFY(!doc.canUndo());  // nothing stored while typing
    QTest::keyClick(cutoff->editor(), Qt::Key_Return);
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("4.75"));
    QCOMPARE(undoSteps(doc), std::size_t{1});
    QCOMPARE(toQString(doc.undoLabel()), QStringLiteral("Set ") + toQString(core::describeField("analyses.contacts.cutoff")));
    typeLine(form.field("output.video.filename")->editor(), "movie.mp4");
    QCOMPARE(undoSteps(doc), std::size_t{2});
    doc.undo();
    QCOMPARE(form.field("output.video.filename")->displayedText(), toQString(text(doc, "output.video.filename")));
    doc.undo();
    QCOMPARE(cutoff->displayedText(), toQString(text(doc, "analyses.contacts.cutoff")));
    QVERIFY(!doc.canUndo());
    QVERIFY(doc.canRedo());
    doc.redo();
    QCOMPARE(cutoff->displayedText(), QStringLiteral("4.75"));
    doc.redo();
    QCOMPARE(form.field("output.video.filename")->displayedText(), QStringLiteral("movie.mp4"));
    QVERIFY(!doc.canRedo());
    // focus leaving a field with typed text stores it once
    auto* title = form.field("visual.title");
    if (form.active) {
        title->editor()->setFocus();
        QTest::keyClicks(title->editor(), "T");
        form.field("visual.center_on")->editor()->setFocus();
        QCOMPARE(text(doc, "visual.title"), std::string("T"));
    }
    // a new project: every widget shows its value
    doc.newProject();
    for (FieldWidget* f : form.fields) {
        if (core::fieldInfo(doc.project(), f->path()) && !qobject_cast<ListEditor*>(f)) {
            QCOMPARE(f->displayedText(), toQString(doc.text(f->path())));
        }
    }
    QVERIFY(!doc.canUndo());
}

void FormsTest::multiLineTextCommitsOnLeaving() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* notes = form.field<MultiLineField>("notes");
    QVERIFY(notes);
    notes->textEdit()->selectAll();
    QTest::keyClicks(notes->textEdit(), "Edited in the window.");
    QTest::keyClick(notes->textEdit(), Qt::Key_Return);
    QTest::keyClicks(notes->textEdit(), "Second line");
    QVERIFY(!doc.canUndo());
    form.context.commitPendingEdits();  // what a save does first
    QCOMPARE(text(doc, "notes"), std::string("Edited in the window.\nSecond line"));
    QCOMPARE(undoSteps(doc), std::size_t{1});
    form.context.commitPendingEdits();
    QCOMPARE(undoSteps(doc), std::size_t{1});
    doc.undo();
    QCOMPARE(notes->displayedText(), QString());
    auto* extra = form.field<MultiLineField>("cleaning.extra_commands");
    QVERIFY(extra);
    if (form.active) {
        extra->textEdit()->setFocus();
        QTest::keyClicks(extra->textEdit(), "image");
        form.field("visual.title")->editor()->setFocus();  // focus out stores
        QCOMPARE(text(doc, "cleaning.extra_commands"), std::string("image"));
    }
}

// The inputs of the GUI harness: Space toggles check boxes, editable combos keep the text typed,
// fixed combos are picked by project value.
void FormsTest::checkBoxesAndEditableCombos() {
    core::ProjectDocument doc;
    Form form(doc);
    auto* pbc = form.field<CheckField>("analyses.use_pbc");
    QVERIFY(pbc);
    const bool before = pbc->checkBox()->isChecked();
    QTest::keyClick(pbc->checkBox(), Qt::Key_Space);
    QCOMPARE(text(doc, "analyses.use_pbc"), std::string(before ? "false" : "true"));
    QCOMPARE(undoSteps(doc), std::size_t{1});
    doc.undo();
    QCOMPARE(pbc->checkBox()->isChecked(), before);
    auto* palette = form.field<EditableComboField>("visual.chain_palette");
    QVERIFY(palette);
    typeLine(palette->editor(), "set2");  // not swapped for the item "Set2"
    QCOMPARE(text(doc, "visual.chain_palette"), std::string("set2"));
    typeLine(palette->editor(), "tab");  // not completed to "tab10"
    QCOMPARE(text(doc, "visual.chain_palette"), std::string("tab"));
    QCOMPARE(palette->comboBox()->count(), 4);  // nothing inserted
    const int dark = palette->comboBox()->findText(QStringLiteral("Dark2"));
    emit palette->comboBox()->activated(dark);
    QCOMPARE(text(doc, "visual.chain_palette"), std::string("Dark2"));
    palette->setChoices({QStringLiteral("a"), QStringLiteral("b")});
    QCOMPARE(palette->displayedText(), QStringLiteral("Dark2"));
    auto* unit = form.field<ComboField>("inputs.time_unit");
    const int us = unit->comboBox()->findData(QStringLiteral("us"));
    QCOMPARE(unit->comboBox()->itemText(us), QStringLiteral("µs"));
    unit->comboBox()->setCurrentIndex(us);
    emit unit->comboBox()->activated(us);
    QCOMPARE(text(doc, "inputs.time_unit"), std::string("us"));
}

void FormsTest::pathFieldsBrowseWithMdatksFilters() {
    TempDir dir;
    core::ProjectDocument doc;
    doc.saveAs(dir / "p.mdatk.toml");
    Form form(doc);
    auto* topology = form.field<PathField>("inputs.topology");
    QVERIFY(topology);
    FileDialogRequest seen;
    const auto previous = setFileDialogFunction([&](QWidget*, const FileDialogRequest& r) {
        seen = r;
        return QStringList{QStringLiteral("/data/run/system.parm7")};
    });
    const QString start = QDir::cleanPath(toQString((dir / doc.text("inputs.topology")).generic_string()));
    topology->browseButton()->click();
    QCOMPARE(seen.title, QStringLiteral("Select topology"));
    QCOMPARE(seen.mode, PathMode::OpenFile);
    QCOMPARE(seen.filters.front(), QStringLiteral("Topology files (*.parm7 *.prmtop *.top *.pdb *.psf *.gro "
                                                  "*.PARM7 *.PRMTOP *.TOP *.PDB *.PSF *.GRO)"));
    QCOMPARE(seen.filters.back(), QStringLiteral("All files (*)"));
    QCOMPARE(seen.start, start);
    QCOMPARE(text(doc, "inputs.topology"), toStdString(QDir::toNativeSeparators("/data/run/system.parm7")));
    QCOMPARE(undoSteps(doc), std::size_t{1});
    // folders
    auto* out = form.field<PathField>("output.out_dir");
    out->browseButton()->click();
    QCOMPARE(seen.mode, PathMode::Folder);
    QVERIFY(seen.filters.isEmpty());
    // typed paths stay as typed (relative)
    typeLine(topology->editor(), "inputs/other.parm7");
    QCOMPARE(text(doc, "inputs.topology"), std::string("inputs/other.parm7"));
    // a list of files: one step for every file picked
    auto* trajectories = form.field<ListEditor>("inputs.trajectories");
    QVERIFY(trajectories && trajectories->browseButton());
    const std::size_t steps = undoSteps(doc);
    const std::size_t count = doc.get("inputs.trajectories").size();
    setFileDialogFunction([&](QWidget*, const FileDialogRequest& r) {
        seen = r;
        return QStringList{QStringLiteral("/a/1.nc"), QStringLiteral("/a/2.nc")};
    });
    trajectories->browseButton()->click();
    QCOMPARE(seen.mode, PathMode::OpenFiles);
    QVERIFY(seen.filters.front().startsWith(QStringLiteral("Trajectory files (*.nc")));
    QCOMPARE(doc.get("inputs.trajectories").size(), count + 2);
    QCOMPARE(undoSteps(doc), steps + 1);
    QCOMPARE(trajectories->rowCount(), static_cast<int>(count + 2));
    auto* external = form.field<ListEditor>("inputs.external");
    QVERIFY(external);
    setFileDialogFunction([&](QWidget*, const FileDialogRequest&) { return QStringList{QStringLiteral("/a/e.dat")}; });
    external->browseButton()->click();
    QCOMPARE(text(doc, "inputs.external.0.path"), toStdString(QDir::toNativeSeparators("/a/e.dat")));
    QCOMPARE(text(doc, "inputs.external.0.x_unit"), std::string("ps"));
    setFileDialogFunction(previous);
}

// Add / remove / move rows (one step each, MDATK's defaults for new rows) and undo of them.
void FormsTest::listEditorRowsAndUndo() {
    core::ProjectDocument doc;
    One one(doc, "groups");
    auto* groups = qobject_cast<ListEditor*>(one.field);
    QVERIFY(groups);
    QCOMPARE(groups->view()->objectName(), QStringLiteral("groups"));
    QCOMPARE(groups->rowCount(), 2);
    QCOMPARE(groups->model()->headerData(3, Qt::Horizontal, ListModel::KeyRole).toString(), QStringLiteral("color"));
    groups->addButton()->click();
    QCOMPARE(groups->rowCount(), 3);
    QCOMPARE(text(doc, "groups.2.name"), std::string("NewGroup"));
    QCOMPARE(text(doc, "groups.2.color"), std::string("#6250d6"));
    QCOMPARE(groups->selectedRow(), 2);
    groups->upButton()->click();
    QCOMPARE(text(doc, "groups.1.name"), std::string("NewGroup"));
    QCOMPARE(groups->selectedRow(), 1);
    groups->selectRow(0);
    groups->removeButton()->click();
    QCOMPARE(text(doc, "groups.0.name"), std::string("NewGroup"));
    QCOMPARE(groups->rowCount(), 2);
    QCOMPARE(undoSteps(doc), std::size_t{3});
    doc.undo();
    QCOMPARE(groups->rowCount(), 3);
    doc.undo();
    doc.undo();
    QCOMPARE(groups->rowCount(), 2);
    QCOMPARE(groups->cellIndex("groups.1.name").data().toString(), QStringLiteral("Target"));
    // the colour cell has a swatch
    QVERIFY(groups->cellIndex("groups.0.color").data(Qt::DecorationRole).isValid());
    // records without a spec take their defaults (RMSD targets, representations)
    One reps(doc, "visual.representations");
    auto* r = qobject_cast<ListEditor*>(reps.field);
    r->addButton()->click();
    QCOMPARE(text(doc, "visual.representations.0.style"), std::string("cartoon"));
    QCOMPARE(r->cellIndex("visual.representations.0.style").data().toString(),
             QStringLiteral("New Cartoon (bases as slabs)"));
    QCOMPARE(r->cellIndex("visual.representations.0.sugars").data(Qt::CheckStateRole).toInt(),
             static_cast<int>(Qt::Unchecked));
    // a cell edit handler (a rename cascade) takes over the edit
    std::vector<std::string> handled;
    groups->setCellEditHandler([&](int row, const std::string& key, const QString& t) {
        handled.push_back(std::to_string(row) + key + toStdString(t));
        doc.renameGroup(doc.text("groups." + std::to_string(row) + ".name"), toStdString(t));
        return true;
    });
    groups->model()->setData(groups->cellIndex("groups.0.name"), QStringLiteral("Apt"), Qt::EditRole);
    QCOMPARE(handled, std::vector<std::string>{"0nameApt"});
    QCOMPARE(text(doc, "groups.0.name"), std::string("Apt"));
    QCOMPARE(groups->cellIndex("groups.0.name").data().toString(), QStringLiteral("Apt"));
}

// AP-15: a cell's typed text is stored when the cell is left by a click elsewhere, not only with
// Enter; Enter stores too; findCell maps cell paths (the harness's keys) to the table.
void FormsTest::tableCellsCommitOnFocusLossAndEnter() {
    core::ProjectDocument doc;
    One one(doc, "groups");
    auto* groups = qobject_cast<ListEditor*>(one.field);
    const CellLocation cell = findCell(one.window.get(), "groups.1.residues");
    QCOMPARE(cell.editor, groups);
    QVERIFY(cell.index.isValid());
    QVERIFY(!findCell(one.window.get(), "groups.7.residues").editor);
    // Enter
    QTableView* view = groups->view();
    view->setCurrentIndex(cell.index);
    view->edit(cell.index);
    auto* editor = qobject_cast<QLineEdit*>(view->indexWidget(cell.index));
    QVERIFY(editor);
    QCOMPARE(editor->objectName(), QStringLiteral("groups.1.residues"));
    editor->selectAll();
    QTest::keyClicks(editor, "A:20-25");
    QTest::keyClick(editor, Qt::Key_Return);
    QTRY_COMPARE(text(doc, "groups.1.residues"), std::string("A:20-25"));  // the view commits queued
    QCOMPARE(undoSteps(doc), std::size_t{1});
    // commitPendingEdits (a save while a cell is open)
    const QModelIndex name = groups->cellIndex("groups.0.name");
    view->edit(name);
    editor = qobject_cast<QLineEdit*>(view->indexWidget(name));
    QVERIFY(editor);
    editor->selectAll();
    QTest::keyClicks(editor, "Apt");
    one.context.commitPendingEdits();
    QCOMPARE(text(doc, "groups.0.name"), std::string("Apt"));
    QCOMPARE(undoSteps(doc), std::size_t{2});
    // a click elsewhere (focus moves out of the cell)
    if (!one.active) {
        QSKIP("the offscreen window could not be activated (focus events need it)");
    }
    const QModelIndex selection = groups->cellIndex("groups.0.selection");
    view->setFocus();
    view->edit(selection);
    editor = qobject_cast<QLineEdit*>(view->indexWidget(selection));
    QVERIFY(editor);
    editor->setFocus();
    QTest::keyClicks(editor, "segid A");
    QTest::mouseClick(one.other, Qt::LeftButton);
    one.other->setFocus();
    QTRY_COMPARE(text(doc, "groups.0.selection"), std::string("segid A"));
    QCOMPARE(undoSteps(doc), std::size_t{3});
    // ... and a click on a button of the editor stores first, then acts
    const QModelIndex colour = groups->cellIndex("groups.1.color");
    view->setFocus();
    view->edit(colour);
    editor = qobject_cast<QLineEdit*>(view->indexWidget(colour));
    editor->setFocus();
    editor->selectAll();
    QTest::keyClicks(editor, "#aa22cc");
    QTest::mouseClick(groups->addButton(), Qt::LeftButton);
    QCOMPARE(text(doc, "groups.1.color"), std::string("#aa22cc"));
    QCOMPARE(groups->rowCount(), 3);
    QCOMPARE(undoSteps(doc), std::size_t{5});
    // a check cell toggles with Space
    One reps(doc, "visual.representations");
    auto* r = qobject_cast<ListEditor*>(reps.field);
    r->addButton()->click();
    const QModelIndex sugars = r->cellIndex("visual.representations.0.sugars");
    r->view()->setFocus();
    r->view()->setCurrentIndex(sugars);
    QTest::keyClick(r->view(), Qt::Key_Space);
    QCOMPARE(text(doc, "visual.representations.0.sugars"), std::string("true"));
    // a combo cell stores the project value of the item picked
    const QModelIndex style = r->cellIndex("visual.representations.0.style");
    r->view()->edit(style);
    auto* combo = qobject_cast<QComboBox*>(r->view()->indexWidget(style));
    QVERIFY(combo);
    combo->setCurrentIndex(combo->findData(QStringLiteral("licorice")));
    reps.context.commitPendingEdits();
    QCOMPARE(text(doc, "visual.representations.0.style"), std::string("licorice"));
}

// MDATK 0.0.4 test_value_widgets_ignore_the_mouse_wheel: every spin box, combo and slider of a
// form ignores the wheel (unless it has the focus: then the wheel is meant for it).
void FormsTest::valueWidgetsIgnoreTheMouseWheel() {
    core::ProjectDocument doc;
    Form form(doc);
    QStringList missing;
    for (QWidget* w : form.window->findChildren<QWidget*>()) {
        const bool valueWidget = qobject_cast<QAbstractSpinBox*>(w) || qobject_cast<QComboBox*>(w) ||
                                 (qobject_cast<QAbstractSlider*>(w) && !qobject_cast<QScrollBar*>(w));
        if (valueWidget && !hasWheelGuard(w)) {
            missing.append(w->objectName());
        }
    }
    QVERIFY2(missing.isEmpty(), qPrintable(missing.join(QStringLiteral(", "))));
    const auto wheel = [](QWidget* target) {
        QWheelEvent event(QPointF(5, 5), target->mapToGlobal(QPointF(5, 5)), QPoint(), QPoint(0, 120), Qt::NoButton,
                          Qt::NoModifier, Qt::NoScrollPhase, false);
        QApplication::sendEvent(target, &event);
    };
    auto* cutoff = form.field<NumberField>("analyses.contacts.cutoff");
    auto* theme = form.field<ComboField>("visual.theme");
    form.field("visual.title")->editor()->setFocus();
    wheel(cutoff->box());
    wheel(theme->comboBox());
    QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("4.5"));
    QCOMPARE(text(doc, "visual.theme"), std::string("light"));
    QVERIFY(!doc.canUndo());
    QCOMPARE(cutoff->box()->focusPolicy(), Qt::StrongFocus);
    if (form.active) {
        cutoff->box()->setFocus();
        QVERIFY(cutoff->box()->hasFocus());
        wheel(cutoff->box());
        QCOMPARE(text(doc, "analyses.contacts.cutoff"), std::string("4.6"));
    }
}

// D10: every editor has an accessible name and takes the keyboard focus.
void FormsTest::everyEditorIsNamedAndReachable() {
    core::ProjectDocument doc;
    Form form(doc);
    for (FieldWidget* f : form.fields) {
        QVERIFY2(!f->editor()->accessibleName().isEmpty(), f->path().c_str());
        QVERIFY2(f->editor()->focusPolicy() & Qt::TabFocus, f->path().c_str());
        QCOMPARE(f->editor()->objectName(), toQString(f->path()));
        QCOMPARE(form.window->findChildren<QWidget*>(toQString(f->path())).size(), qsizetype{1});
    }
    for (QAbstractButton* b : form.window->findChildren<QAbstractButton*>()) {
        if (!b->objectName().contains(QLatin1Char(':'))) {
            continue;  // Qt's own (a line edit's icon)
        }
        QVERIFY2(!b->accessibleName().isEmpty() || !b->text().isEmpty(), qPrintable(b->objectName()));
        QVERIFY2(b->focusPolicy() & Qt::TabFocus, qPrintable(b->objectName()));
    }
    QCOMPARE(form.field("analyses.contacts.cutoff")->label(), QStringLiteral("Cutoff"));
}

QTEST_MAIN(FormsTest)
#include "test_forms.moc"
