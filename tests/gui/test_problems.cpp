// The problems presentation (src/app/gui/problems, milestone 5 W1-D), Qt Test on the offscreen
// platform: every problem of every validate case (tests/data/project-model/validate_cases.json)
// lands on a field of a form built from the field factory and the panel shows it with its
// location; page badges count each page's problems; a click on a row jumps to the field; typed
// text is checked as you type, once per burst (MDATK 0.1.1's CHECK_DELAY_MS); the field widgets'
// own warnings (the range hint) are listed as warnings.

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <QApplication>
#include <QLineEdit>
#include <QListWidget>
#include <QScrollArea>
#include <QSignalSpy>
#include <QStackedWidget>
#include <QTableView>
#include <QVBoxLayout>
#include <QtTest>

#include "gui/display_path.hpp"
#include "gui/forms/field_factory.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/problems/problems_panel.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/project_io.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::app::forms;
using namespace mdbridge::app::problems;
using mdbridge::test::TempDir;
namespace fs = std::filesystem;

namespace {

const fs::path kData = MDBRIDGE_TEST_DATA;

class ScopedCwd {
public:
    explicit ScopedCwd(const fs::path& path) : m_old(fs::current_path()) { fs::current_path(path); }
    ~ScopedCwd() { fs::current_path(m_old); }
    ScopedCwd(const ScopedCwd&) = delete;
    ScopedCwd& operator=(const ScopedCwd&) = delete;

private:
    fs::path m_old;
};

std::vector<std::string> strings(const core::Value& value) {
    std::vector<std::string> out;
    for (const core::Value& v : *value.get<core::Array>()) {
        out.push_back(v.stringOr());
    }
    return out;
}

// A test-only window as W1-A's will be: one page per side-bar page (a scroll area of the factory's
// field widgets of that page), the problems panel and a badge per page; the panel's reveal function
// shows the field's page and reveals the field there (what the window does with a page's
// revealField).
struct Host {
    explicit Host(core::ProjectDocument& doc) : context(doc) {
        window = std::make_unique<QWidget>();
        auto* outer = new QHBoxLayout(window.get());
        stack = new QStackedWidget(window.get());
        auto* side = new QVBoxLayout;
        for (const core::ProblemPage& page : core::kProblemPages) {
            const QString id = toQString(page.id);
            auto* scroll = new QScrollArea;
            scroll->setObjectName(QStringLiteral("page:") + id);
            scroll->setWidgetResizable(true);
            auto* body = new QWidget;
            layouts[id] = new QVBoxLayout(body);
            scroll->setWidget(body);
            stack->addWidget(scroll);
            pages[id] = scroll;
            auto* badge = new ProblemBadge(context, id, window.get());
            badges[id] = badge;
            side->addWidget(badge);
        }
        for (const std::string& path : formFieldPaths(doc.project())) {
            const QString page = toQString(core::problemPage(path));
            if (!layouts.count(page)) {
                continue;  // no page shows it (the planner's workflow)
            }
            if (FieldWidget* field = createFieldWidget(context, path, layouts[page]->parentWidget())) {
                layouts[page]->addWidget(field);
            }
        }
        panel = new ProblemsPanel(context, window.get());
        panel->setRevealFunction([this](const QString& field) {
            const QString page = toQString(core::problemPage(toStdString(field)));
            if (!pages.count(page)) {
                return false;
            }
            stack->setCurrentWidget(pages[page]);
            return revealFieldIn(pages[page], field);
        });
        side->addWidget(panel, 1);
        outer->addLayout(side);
        outer->addWidget(stack, 1);
        window->resize(1000, 700);
        window->show();
        window->activateWindow();
        active = QTest::qWaitForWindowActive(window.get());
        context.refreshProblems();
        settle();
    }

    // Lets the coalesced refreshes of the panel and badges run.
    void settle() {
        QCoreApplication::processEvents();
        panel->refresh();
        for (auto& [id, badge] : badges) {
            badge->refresh();
        }
    }

    const core::Problem* shown(const std::string& field, const std::string& message) const {
        for (const core::Problem& p : panel->problems()) {
            if (p.field == field && p.message == message) {
                return &p;
            }
        }
        return nullptr;
    }

    QListWidgetItem* row(const std::string& field, const std::string& message) const {
        QListWidget* list = panel->list();
        for (int i = 0; i < list->count(); ++i) {
            const core::Problem& p = panel->problems()[static_cast<std::size_t>(i)];
            if (p.field == field && p.message == message) {
                return list->item(i);
            }
        }
        return nullptr;
    }

    FormContext context;
    std::unique_ptr<QWidget> window;
    QStackedWidget* stack = nullptr;
    ProblemsPanel* panel = nullptr;
    std::map<QString, QScrollArea*> pages;
    std::map<QString, QVBoxLayout*> layouts;
    std::map<QString, ProblemBadge*> badges;
    bool active = false;
};

// The path of the field widget or list editor the focus is in ("" when none).
std::string focusedPath() {
    for (QWidget* w = QApplication::focusWidget(); w; w = w->parentWidget()) {
        if (auto* field = qobject_cast<FieldWidget*>(w)) {
            return field->path();
        }
    }
    return {};
}

// The field path the focus shows: a field widget's, or a list editor's current cell.
std::string focusedField() {
    for (QWidget* w = QApplication::focusWidget(); w; w = w->parentWidget()) {
        if (auto* list = qobject_cast<ListEditor*>(w)) {
            const QModelIndex index = list->view()->currentIndex();
            return index.isValid() ? list->model()->cellPath(index.row(), index.column()) : list->path();
        }
        if (auto* field = qobject_cast<FieldWidget*>(w)) {
            return field->path();
        }
    }
    return {};
}

// The line edit of a one-line field (a line edit, spin box or editable combo).
QLineEdit* lineOf(FieldWidget* field) {
    if (!field) {
        return nullptr;
    }
    if (auto* edit = qobject_cast<QLineEdit*>(field->editor())) {
        return edit;
    }
    return field->editor()->findChild<QLineEdit*>();
}

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

bool exactWidget(QWidget* root, const std::string& field) {
    return root->findChild<FieldWidget*>(toQString(field) + QStringLiteral(":field")) ||
           findCell(root, field).editor;
}

}  // namespace

class ProblemsTest : public QObject {
    Q_OBJECT

private slots:
    void everyValidateProblemLandsOnItsField();
    void badgesCountEachPage();
    void clickingARowJumpsToTheField();
    void typedTextIsCheckedOncePerBurst();
    void escapeDropsTheTypedProblem();
    void aDragIsCheckedWhenItPauses();
    void rangeHintWarningsAreListed();
    void summariesAndCounts();
};

// Every validate() problem of every case: shown by the panel with describeField's location, on an
// existing field (a field widget or list cell of exactly that path), and revealing it focuses that
// field on its page.
void ProblemsTest::everyValidateProblemLandsOnItsField() {
    const core::Value vectors = mdbridge::test::readJson(kData / "project-model" / "validate_cases.json");
    TempDir dir;
    for (const std::string& name : strings(*vectors.find("files"))) {
        mdbridge::test::writeFile(dir / name, "");
    }
    for (const std::string& name : strings(*vectors.find("folders"))) {
        fs::create_directory(dir / name);
    }
    const ScopedCwd cwd(dir.path());
    std::size_t landed = 0;
    std::size_t cases = 0;
    for (const core::Value& c : *vectors.find("cases")->get<core::Array>()) {
        const std::string name = c.find("name")->stringOr();
        core::ProjectDocument doc;
        doc.replace(core::projectFromPlain(mdbridge::test::tagged(*c.find("input"))));
        Host host(doc);
        QVERIFY2(host.active, name.c_str());
        ++cases;
        const std::vector<core::Problem> expected = doc.problems();
        for (const core::Problem& p : expected) {
            const std::string what = name + ": " + p.field + ": " + p.message;
            const core::Problem* shown = host.shown(p.field, p.message);
            QVERIFY2(shown, what.c_str());
            QCOMPARE(shown->level, p.level);
            QCOMPARE(shown->location, core::describeField(p.field, &doc.project()));
            QListWidgetItem* item = host.row(p.field, p.message);
            QVERIFY2(item, what.c_str());
            // the row shows the message's paths natively (backslashes on Windows, WN-12)
            QCOMPARE(item->text(),
                     toQString(p.location) + QLatin1Char('\n') + app::gui::displayPathsIn(toQString(p.message)));
            QCOMPARE(item->data(ProblemsPanel::FieldRole).toString(), toQString(p.field));
            if (p.level != "error") {
                continue;  // the layout warnings are at "layout.panels", the list editor of the slots
            }
            const QString page = toQString(p.page());
            QVERIFY2(host.pages.count(page), what.c_str());
            QVERIFY2(exactWidget(host.pages[page], p.field), what.c_str());
            QVERIFY2(host.panel->revealProblem(toQString(p.field)), what.c_str());
            QCOMPARE(host.stack->currentWidget(), host.pages[page]);
            QVERIFY2(focusedField() == p.field || focusedPath() == p.field, (what + " -> " + focusedField()).c_str());
            ++landed;
        }
        // the panel shows exactly the document's problems (no field widget has one of its own here
        // that the checks do not have, except the range hint)
        for (const core::Problem& p : host.panel->problems()) {
            const bool known = std::any_of(expected.begin(), expected.end(), [&](const core::Problem& e) {
                return e.field == p.field && e.message == p.message;
            });
            QVERIFY2(known || p.level == "warning", (name + ": " + p.field + ": " + p.message).c_str());
        }
    }
    QCOMPARE(cases, std::size_t{120});
    QVERIFY(landed > 150);
}

// The badges count the problems of their page (errors red, warnings only amber, hidden at none).
void ProblemsTest::badgesCountEachPage() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;  // the defaults: no topology, no trajectory
    Host host(doc);
    QVERIFY(host.active);
    QCOMPARE(host.badges["input"]->count().errors, 2);
    QCOMPARE(host.badges["input"]->text(), QStringLiteral("2"));
    QVERIFY(host.badges["input"]->isVisible());
    QCOMPARE(host.badges["input"]->property("level").toString(), QStringLiteral("error"));
    QVERIFY(host.badges["input"]->toolTip().contains(QStringLiteral("No topology file selected.")));
    for (const char* page : {"groups", "analyses", "structure", "layout", "outputs"}) {
        QCOMPARE(host.badges[page]->count().total(), host.panel->pageCount(page).total());
    }
    // a group without residues: one more on Residue groups
    QSignalSpy spy(host.badges["groups"], &ProblemBadge::countChanged);
    doc.set("groups.0.residues", core::Value(std::string()));
    host.settle();
    const ProblemCount groups = countOn(host.context.allProblems(), "groups");
    QVERIFY(groups.total() > 0);
    QCOMPARE(host.badges["groups"]->count().total(), groups.total());
    QCOMPARE(host.badges["groups"]->text(), QString::number(groups.total()));
    QCOMPARE(spy.count(), 1);
    // the badges add up to the list
    int sum = 0;
    for (auto& [id, badge] : host.badges) {
        sum += badge->count().total();
    }
    QCOMPARE(sum, host.panel->count().total());
    // a fixed project: no badge
    doc.undo();
    doc.set("inputs.topology", core::Value(std::string("top.parm7")));
    host.settle();
    mdbridge::test::writeFile(dir / "top.parm7", "");
    mdbridge::test::writeFile(dir / "a.nc", "");
    doc.set("inputs.trajectories", core::Value(core::Array{core::Value(std::string("a.nc"))}));
    host.settle();
    QCOMPARE(host.badges["input"]->count().total(), 0);
    QVERIFY(!host.badges["input"]->isVisible());
    // clicking a badged page cycles through its problems (MDATK 0.1.1)
    doc.set("inputs.step", core::Value(std::int64_t{0}));
    doc.set("inputs.start", core::Value(std::int64_t{-1}));
    host.settle();
    const int n = host.panel->pageCount("input").total();
    QVERIFY(n >= 2);
    const core::Problem* first = host.panel->nextOnPage("input");
    const core::Problem* second = host.panel->nextOnPage("input");
    QVERIFY(first && second && first != second);
    for (int k = 2; k < n; ++k) {
        host.panel->nextOnPage("input");
    }
    QCOMPARE(host.panel->nextOnPage("input"), first);
    QCOMPARE(host.panel->nextOnPage("groups"), nullptr);
}

void ProblemsTest::clickingARowJumpsToTheField() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;
    Host host(doc);
    QVERIFY(host.active);
    host.stack->setCurrentWidget(host.pages["outputs"]);
    QSignalSpy spy(host.panel, &ProblemsPanel::problemActivated);
    QListWidgetItem* item = host.row("inputs.trajectories", "No trajectory file selected.");
    QVERIFY(item);
    QListWidget* list = host.panel->list();
    list->scrollToItem(item);
    QTest::mouseClick(list->viewport(), Qt::LeftButton, {}, list->visualItemRect(item).center());
    QCOMPARE(spy.count(), 1);
    QCOMPARE(spy.front().at(0).toString(), QStringLiteral("inputs.trajectories"));
    QCOMPARE(spy.front().at(1).toString(), QStringLiteral("input"));
    QCOMPARE(host.stack->currentWidget(), host.pages["input"]);
    QCOMPARE(focusedPath(), std::string("inputs.trajectories"));
    // the field is in view
    QWidget* field = host.pages["input"]->findChild<FieldWidget*>(QStringLiteral("inputs.trajectories:field"));
    QVERIFY(field);
    const QRect inView(field->mapTo(host.pages["input"]->viewport(), QPoint(0, 0)), field->size());
    QVERIFY(host.pages["input"]->viewport()->rect().intersects(inView));
    // a line edit's text is not selected by the jump
    doc.set("inputs.topology", core::Value(std::string("missing.parm7")));
    host.settle();
    QListWidgetItem* topology = nullptr;
    for (int i = 0; i < list->count(); ++i) {
        if (list->item(i)->data(ProblemsPanel::FieldRole).toString() == QStringLiteral("inputs.topology")) {
            topology = list->item(i);
        }
    }
    QVERIFY(topology);
    QTest::keyClick(list, Qt::Key_Home);
    list->setCurrentItem(topology);
    QTest::keyClick(list, Qt::Key_Return);
    auto* edit = qobject_cast<QLineEdit*>(QApplication::focusWidget());
    QVERIFY(edit);
    QCOMPARE(edit->objectName(), QStringLiteral("inputs.topology"));
    QVERIFY(!edit->hasSelectedText());
    // without a reveal function nothing happens
    host.panel->setRevealFunction({});
    QVERIFY(!host.panel->revealProblem(QStringLiteral("inputs.topology")));
}

// Typing is checked kCheckDelayMs after the last key: one check for a burst, on the text as typed
// (not yet stored), and the field is flagged with it.
void ProblemsTest::typedTextIsCheckedOncePerBurst() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;
    Host host(doc);
    QVERIFY(host.active);
    auto* stop = host.pages["input"]->findChild<FieldWidget*>(QStringLiteral("inputs.stop:field"));
    QLineEdit* stopEdit = lineOf(stop);
    QVERIFY(stopEdit);
    stop->setFocus();
    QSignalSpy spy(host.panel, &ProblemsPanel::problemsShown);
    for (const char key : std::string("abc")) {
        QTest::keyClicks(stopEdit, QString(QLatin1Char(key)));
        QVERIFY(host.context.checkPending());
        QTest::qWait(FormContext::kCheckDelayMs / 4);
    }
    QCOMPARE(spy.count(), 0);  // not yet: still typing
    QVERIFY(!host.shown("inputs.stop", "'abc' is not a whole number"));
    QTRY_VERIFY_WITH_TIMEOUT(host.shown("inputs.stop", "'abc' is not a whole number"), 3000);
    QCOMPARE(spy.count(), 1);
    QVERIFY(!doc.isModified());  // checked, not stored
    QCOMPARE(stop->flag().first, QStringLiteral("'abc' is not a whole number"));
    const core::Problem* typed = host.shown("inputs.stop", "'abc' is not a whole number");
    QCOMPARE(typed->location, std::string("Input & cleaning › Frames › Stop"));
    QCOMPARE(typed->level, std::string("error"));
    // a typed value the checks refuse: Step 0
    QLineEdit* step = lineOf(host.pages["input"]->findChild<FieldWidget*>(QStringLiteral("inputs.step:field")));
    QVERIFY(step);
    step->selectAll();
    QTest::keyClicks(step, QStringLiteral("0"));
    QTRY_VERIFY_WITH_TIMEOUT(countOn(host.panel->problems(), "input").errors >= 4, 3000);
    QVERIFY(std::any_of(host.panel->problems().begin(), host.panel->problems().end(),
                        [](const core::Problem& p) { return p.field == "inputs.step"; }));
    QVERIFY(!doc.isModified());
}

void ProblemsTest::escapeDropsTheTypedProblem() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;
    Host host(doc);
    QVERIFY(host.active);
    auto* stop = host.pages["input"]->findChild<FieldWidget*>(QStringLiteral("inputs.stop:field"));
    QLineEdit* edit = lineOf(stop);
    QVERIFY(edit);
    stop->setFocus();
    QTest::keyClicks(edit, QStringLiteral("x"));
    QTRY_VERIFY_WITH_TIMEOUT(host.shown("inputs.stop", "'x' is not a whole number"), 3000);
    QTest::keyClick(edit, Qt::Key_Escape);
    QTRY_VERIFY_WITH_TIMEOUT(!host.shown("inputs.stop", "'x' is not a whole number"), 3000);
    QCOMPARE(stop->flag().first, QString());
}

// The changes of an undo group (a drag: one per mouse move) are not checked one by one: the last
// check stands until the drag pauses for kCheckDelayMs.
void ProblemsTest::aDragIsCheckedWhenItPauses() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;
    Host host(doc);
    QVERIFY(host.active);
    const auto stepProblems = [&] { return host.context.problemsAt("inputs.step").size(); };
    QCOMPARE(stepProblems(), std::size_t{0});
    doc.beginGroup("Drag");
    for (std::int64_t v : {3, 2, 1, 0}) {
        doc.set("inputs.step", core::Value(v));
        QCoreApplication::processEvents();
        QCOMPARE(stepProblems(), std::size_t{0});  // not checked per move
        QVERIFY(host.context.checkPending());
    }
    QTRY_VERIFY_WITH_TIMEOUT(stepProblems() > 0, 3000);  // checked once it paused
    doc.endGroup();
    QCOMPARE(undoSteps(doc), std::size_t{1});
    // outside a group every change is checked at once
    doc.undo();
    QCOMPARE(stepProblems(), std::size_t{0});
}

// The range hint's warning (a field widget's own problem, which validate() allows) is listed as a
// warning at the field, counted by its page's (amber) badge, and goes away with the value.
void ProblemsTest::rangeHintWarningsAreListed() {
    TempDir dir;
    const ScopedCwd cwd(dir.path());
    core::ProjectDocument doc;
    Host host(doc);
    QVERIFY(host.active);
    QCOMPARE(host.badges["analyses"]->count().total(), 0);
    doc.set("analyses.contacts.cutoff", core::Value(25.0));
    host.settle();
    QVERIFY(host.context.problemsAt("analyses.contacts.cutoff").empty());  // validate() allows it
    const core::Problem* warning = nullptr;
    for (const core::Problem& p : host.panel->problems()) {
        if (p.field == "analyses.contacts.cutoff") {
            warning = &p;
        }
    }
    QVERIFY(warning);
    QCOMPARE(warning->level, std::string("warning"));
    QVERIFY(QString::fromStdString(warning->message).contains(QStringLiteral("outside the usual range")));
    QCOMPARE(warning->location, core::describeField("analyses.contacts.cutoff", &doc.project()));
    QCOMPARE(host.badges["analyses"]->count().warnings, 1);
    QCOMPARE(host.badges["analyses"]->property("level").toString(), QStringLiteral("warning"));
    QVERIFY(host.panel->revealProblem(QStringLiteral("analyses.contacts.cutoff")));
    QCOMPARE(focusedPath(), std::string("analyses.contacts.cutoff"));
    // while typing, the hint follows the text
    auto* cutoff = qobject_cast<NumberField*>(QApplication::focusWidget()->parentWidget());
    QVERIFY(cutoff);
    cutoff->box()->edit()->selectAll();
    QTest::keyClicks(cutoff->box()->edit(), QStringLiteral("4.5"));
    host.settle();
    QCOMPARE(host.badges["analyses"]->count().total(), 0);
    QTest::keyClick(cutoff->box()->edit(), Qt::Key_Return);
    QCOMPARE(doc.text("analyses.contacts.cutoff"), std::string("4.5"));
    host.settle();
    QCOMPARE(host.badges["analyses"]->count().total(), 0);
    QVERIFY(!host.badges["analyses"]->isVisible());
}

void ProblemsTest::summariesAndCounts() {
    QCOMPARE(fieldParents("groups.1.residues"),
             (std::vector<std::string>{"groups.1.residues", "groups.1", "groups"}));
    QVERIFY(fieldParents("").empty());
    QCOMPARE(countText({2, 1}), QStringLiteral("2 errors, 1 warning"));
    QCOMPARE(countText({1, 0}), QStringLiteral("1 error"));
    QCOMPARE(countText({0, 3}), QStringLiteral("3 warnings"));
    std::vector<core::Problem> problems;
    for (int i = 0; i < 8; ++i) {
        problems.push_back(core::locatedProblem("inputs.stop", "bad " + std::to_string(i), nullptr));
    }
    const QString summary = summaryText(problems);
    QVERIFY(summary.startsWith(QStringLiteral("• Frames › Stop: bad 0\n")));
    QVERIFY(summary.endsWith(QStringLiteral("\n… and 2 more")));
    QCOMPARE(summary.count(QLatin1Char('\n')), 6);
    QCOMPARE(countOn(problems, "input").errors, 8);
    QCOMPARE(countOn(problems, "groups").total(), 0);
}

QTEST_MAIN(ProblemsTest)
#include "test_problems.moc"
