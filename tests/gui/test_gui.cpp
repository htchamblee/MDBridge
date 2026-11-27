// mdbridge_gui_tests: the main window driven by action scripts (tests/data/gui, format in
// tests/data/README.md "gui"; tests/parity/make_gui_harness_cases.py and the per-page generators),
// the saved files compared with MDATK's, screenshots for CI artifacts. Qt Test on the "offscreen"
// platform (Linux and Windows, no display), the software 3D view.

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFile>
#include <QLineEdit>
#include <QLocale>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSignalSpy>
#include <QSpinBox>
#include <QTemporaryDir>
#include <QTest>
#include <QVBoxLayout>

#include "gui/forms/field_factory.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/forms/list_editor.hpp"
#include "gui/main_window.hpp"
#include "gui/preview_widget.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"
#include "gui/jobs/log_view.hpp"
#include "mdbridge/core/project_document.hpp"
#include "mdbridge/core/settings.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;
using guitest::CaseError;

namespace {

const fs::path kData = core::pathFromUtf8(MDBRIDGE_TEST_DATA);

fs::path screenshots() { return guitest::screenshotDir(core::pathFromUtf8(MDBRIDGE_GUI_SCREENSHOT_DIR)); }

template <class F>
std::string errorOf(F&& f) {
    try {
        f();
    } catch (const CaseError& e) {
        return e.what();
    }
    return {};
}

}  // namespace

class GuiTests : public QObject {
    Q_OBJECT

private slots:
    void platform();
    void findsAndSetsFields();
    void entersFormFieldsAndTableCells();
    void comparesTexts();
    void readsCaseFiles();
    void cases_data();
    void cases();
};

void GuiTests::platform() {
    // QT_QPA_PLATFORM=offscreen unless the caller chose another platform (main below)
    if (qEnvironmentVariable("MDBRIDGE_GUI_TEST_PLATFORM").isEmpty()) {
        QCOMPARE(QGuiApplication::platformName(), QStringLiteral("offscreen"));
    }
}

// The field lookup and entry on stand-in widgets named like the settings pages' fields will be.
void GuiTests::findsAndSetsFields() {
    QWidget root;
    auto* layout = new QVBoxLayout(&root);
    auto* text = new QLineEdit(QStringLiteral("old"));
    text->setObjectName(QStringLiteral("workflow.project_name"));
    auto* real = new QDoubleSpinBox;
    real->setObjectName(QStringLiteral("analyses.contacts.cutoff"));
    real->setRange(0, 1000);
    real->setDecimals(3);
    auto* whole = new QSpinBox;
    whole->setObjectName(QStringLiteral("inputs.step"));
    whole->setRange(1, 1000);
    auto* check = new QCheckBox;
    check->setObjectName(QStringLiteral("video.enabled"));
    check->setChecked(true);
    auto* combo = new QComboBox;
    combo->setObjectName(QStringLiteral("inputs.time_unit"));
    for (const char* unit : {"ps", "ns", "us"}) {  // label and project value, as the pages' drop-downs
        combo->addItem(QStringLiteral("%1 (label)").arg(QLatin1String(unit)), QString::fromLatin1(unit));
    }
    auto* notes = new QPlainTextEdit(QStringLiteral("old notes"));
    notes->setObjectName(QStringLiteral("notes"));
    auto* button = new QPushButton(QStringLiteral("Add"));
    button->setObjectName(QStringLiteral("groups.add"));
    QSignalSpy clicked(button, &QPushButton::clicked);
    for (QWidget* w : {static_cast<QWidget*>(text), static_cast<QWidget*>(real), static_cast<QWidget*>(whole),
                       static_cast<QWidget*>(check), static_cast<QWidget*>(combo), static_cast<QWidget*>(notes),
                       static_cast<QWidget*>(button)}) {
        layout->addWidget(w);
    }
    root.show();

    guitest::setField(guitest::findField(&root, QStringLiteral("workflow.project_name")), "apt1");
    QCOMPARE(text->text(), QStringLiteral("apt1"));
    guitest::setField(guitest::findField(&root, QStringLiteral("analyses.contacts.cutoff")), 4.25);
    QCOMPARE(real->value(), 4.25);
    guitest::setField(guitest::findField(&root, QStringLiteral("inputs.step")), 5);
    QCOMPARE(whole->value(), 5);
    guitest::setField(guitest::findField(&root, QStringLiteral("video.enabled")), false);
    QVERIFY(!check->isChecked());
    guitest::setField(guitest::findField(&root, QStringLiteral("inputs.time_unit")), "ns");
    QCOMPARE(combo->currentData().toString(), QStringLiteral("ns"));
    guitest::setField(guitest::findField(&root, QStringLiteral("notes")), "line 1\nline 2");
    QCOMPARE(notes->toPlainText(), QStringLiteral("line 1\nline 2"));

    QTest::mouseClick(qobject_cast<QPushButton*>(guitest::findField(&root, QStringLiteral("groups.add"))), Qt::LeftButton);
    QCOMPARE(clicked.count(), 1);

    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::findField(&root, QStringLiteral("missing.key")); })),
             QStringLiteral("no widget is named 'missing.key'"));
    auto* twin = new QLineEdit;
    twin->setObjectName(QStringLiteral("notes"));
    layout->addWidget(twin);
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::findField(&root, QStringLiteral("notes")); })),
             QStringLiteral("2 widgets are named 'notes'"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::setField(combo, "fs"); })),
             QStringLiteral("inputs.time_unit: the drop-down has no item whose value is 'fs'"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::setField(combo, "ns (label)"); })),
             QStringLiteral("inputs.time_unit: the drop-down has no item whose value is 'ns (label)'"));
    // "type": "" clears the field
    guitest::typeField(text, QString());
    QCOMPARE(text->text(), QString());
    guitest::typeField(notes, QString());
    QCOMPARE(notes->toPlainText(), QString());
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::setField(check, 1); })),
             QStringLiteral("video.enabled: a check box takes true or false"));
}

// The harness on the real field widgets (W0-C) bound to a document (W0-B): "set" picks drop-down
// items by project value, "type" with "" clears, table cells by field path, check cells with Space,
// the list editors' rows by "table".
void GuiTests::entersFormFieldsAndTableCells() {
    core::ProjectDocument doc;  // the defaults: groups Target and Aptamer, two trajectories
    app::forms::FormContext context(doc);
    QWidget root;
    auto* layout = new QVBoxLayout(&root);
    for (const char* path : {"notes", "inputs.stop", "inputs.dt_ps", "inputs.time_unit", "visual.center_on", "groups",
                             "inputs.trajectories", "visual.representations"}) {
        app::forms::FieldWidget* field = app::forms::createFieldWidget(context, path, &root);
        QVERIFY2(field, path);
        layout->addWidget(field);
    }
    root.resize(700, 900);
    root.show();
    QVERIFY(QTest::qWaitForWindowExposed(&root));
    const auto text = [&](const char* path) { return QString::fromStdString(doc.text(path)); };

    guitest::setKey(&root, QStringLiteral("inputs.time_unit"), "ns");
    QCOMPARE(text("inputs.time_unit"), QStringLiteral("ns"));
    guitest::typeKey(&root, QStringLiteral("inputs.stop"), QStringLiteral("250"));
    QCOMPARE(text("inputs.stop"), QStringLiteral("250"));
    guitest::typeKey(&root, QStringLiteral("inputs.stop"), QString());  // "" clears: unset
    QCOMPARE(text("inputs.stop"), QString());
    QVERIFY(doc.get("inputs.stop").isNull());
    guitest::typeKey(&root, QStringLiteral("inputs.dt_ps"), QStringLiteral("0.0123456789"));
    QCOMPARE(text("inputs.dt_ps"), QStringLiteral("0.0123456789"));
    guitest::typeKey(&root, QStringLiteral("notes"), QStringLiteral("Edited.\nSecond line"));
    QCOMPARE(text("notes"), QStringLiteral("Edited.\nSecond line"));
    guitest::typeKey(&root, QStringLiteral("visual.center_on"), QStringLiteral("B:12"));
    QCOMPARE(text("visual.center_on"), QStringLiteral("B:12"));

    // rows and cells (tests/data/gui/pages/groups_page.json's steps)
    guitest::tableAction(&root, QStringLiteral("groups"), "add", -1, 0);
    QCOMPARE(doc.get("groups").size(), std::size_t{3});
    guitest::setKey(&root, QStringLiteral("groups.2.name"), "Loop");
    guitest::setKey(&root, QStringLiteral("groups.2.residues"), "A:20-25");
    guitest::typeKey(&root, QStringLiteral("groups.2.color"), QStringLiteral("#aa22cc"));
    QCOMPARE(text("groups.2.name"), QStringLiteral("Loop"));
    QCOMPARE(text("groups.2.residues"), QStringLiteral("A:20-25"));
    QCOMPARE(text("groups.2.color"), QStringLiteral("#aa22cc"));
    guitest::typeKey(&root, QStringLiteral("groups.2.selection"), QString());
    QCOMPARE(text("groups.2.selection"), QString());
    guitest::tableAction(&root, QStringLiteral("groups"), "move", 2, -1);
    QCOMPARE(text("groups.1.name"), QStringLiteral("Loop"));
    guitest::tableAction(&root, QStringLiteral("groups"), "remove", 0, 0);
    QCOMPARE(doc.get("groups").size(), std::size_t{2});
    QCOMPARE(text("groups.0.name"), QStringLiteral("Loop"));
    const std::size_t trajectories = doc.get("inputs.trajectories").size();
    guitest::tableAction(&root, QStringLiteral("inputs.trajectories"), "add", -1, 0);
    const QString last = QStringLiteral("inputs.trajectories.%1").arg(trajectories);
    guitest::setKey(&root, last, "md3.nc");
    QCOMPARE(QString::fromStdString(doc.text(last.toStdString())), QStringLiteral("md3.nc"));
    guitest::tableAction(&root, QStringLiteral("visual.representations"), "add", -1, 0);
    guitest::setKey(&root, QStringLiteral("visual.representations.0.style"), "licorice");  // a drop-down cell
    QCOMPARE(text("visual.representations.0.style"), QStringLiteral("licorice"));
    guitest::setKey(&root, QStringLiteral("visual.representations.0.sugars"), true);  // a check cell
    QCOMPARE(text("visual.representations.0.sugars"), QStringLiteral("true"));
    guitest::setKey(&root, QStringLiteral("visual.representations.0.sugars"), false);
    QCOMPARE(text("visual.representations.0.sugars"), QStringLiteral("false"));

    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::setKey(&root, QStringLiteral("groups.9.name"), "x"); })),
             QStringLiteral("no widget or table cell is named 'groups.9.name'"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::tableAction(&root, QStringLiteral("groups"), "remove", 5, 0); })),
             QStringLiteral("groups: there is no row 5 (2 rows)"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::tableAction(&root, QStringLiteral("notes"), "add", -1, 0); })),
             QStringLiteral("no list editor edits 'notes'"));
    QCOMPARE(QString::fromStdString(
                 errorOf([&] { guitest::setKey(&root, QStringLiteral("visual.representations.0.style"), "glow"); })),
             QStringLiteral("visual.representations.0.style: the drop-down has no item whose value is 'glow'"));
}

void GuiTests::comparesTexts() {
    QCOMPARE(guitest::compareText("a\nb\n", "a\nb\n"), std::string());
    QCOMPARE(guitest::compareText("a\nc\n", "a\nb\n"), std::string("line 2: expected 'b', got 'c'"));
    QCOMPARE(guitest::compareText("a\n", "a\nb\n"), std::string("line 2: expected 'b', got the end of the file"));
    QCOMPARE(guitest::compareText("a", "a\n"), std::string("the texts differ in their final line end"));
    // a fitted case's texts (D10): decimals within their last printed digit or the file's
    // tolerance, the spaces before them free (a sign moves a PDB column), whole numbers exact
    QVERIFY(guitest::closeText("CC 0.912 at step 3", "CC 0.913 at step 3"));
    QVERIFY(!guitest::closeText("CC 0.912 at step 3", "CC 0.914 at step 3"));
    QVERIFY(!guitest::closeText("CC 0.912 at step 3", "CC 0.912 at step 4"));
    QVERIFY(guitest::closeText("ATOM      1  P     G A   1      -0.001  12.500", "ATOM      1  P     G A   1       0.004  12.500",
                               0.01));
    QVERIFY(!guitest::closeText("x   1.000", "x   1.020", 0.01));
    QVERIFY(guitest::closeText(" 1.5", "-1.5", 3.0));
    QCOMPARE(guitest::compareFittedText("a 1.25\nb\n", "a 1.24\nb\n"), std::string());
    QCOMPARE(guitest::compareFittedText("a 1.25\nb", "a 1.24\nb\n"), std::string("the texts differ in their final line end"));
    QCOMPARE(guitest::compareFittedText("a 1.27\n", "a 1.24\n"),
             std::string("line 1: expected 'a 1.24', got 'a 1.27' (decimals compared within their last digit or 0.0)"));
    // a project file: the [density] table by D10's angle / shift / CC tolerances, the rest by line
    const std::string want = "[density]\ntransform = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0]\n"
                             "anchor = [[0.0, 0.0, 0.0]]\nfit_cc = 0.9\nfit_info = \"CC 0.900\"\n";
    std::string got = want;
    got.replace(got.find("1.0, 2.0"), 8, "1.1, 2.0");
    got.replace(got.find("0.9\n"), 4, "0.9005\n");
    got.replace(got.find("CC 0.900"), 8, "CC 0.901");
    QCOMPARE(guitest::compareFittedFile("p.mdatk.toml", got, want, 0.0), std::string());
    got.replace(got.find("1.1, 2.0"), 8, "1.3, 2.0");
    QVERIFY(guitest::compareFittedFile("p.mdatk.toml", got, want, 0.0).rfind("the [density] table differs", 0) == 0);
    // Windows: the folder with either separator, the path after it with '/' up to where it ends
    QCOMPARE(guitest::replaceFolder("Saved project C:/work/tmp\\case\\sub\\copy.mdatk.toml", "C:/work/tmp/case",
                                    "{case_dir}", true),
             std::string("Saved project {case_dir}/sub/copy.mdatk.toml"));
    QCOMPARE(guitest::replaceFolder("ERROR: Could not open C:\\work\\case\\a.toml: C:\\work\\case\\a.toml: bad '\\n'",
                                    "C:/work/case", "{case_dir}", true),
             std::string("ERROR: Could not open {case_dir}/a.toml: {case_dir}/a.toml: bad '\\n'"));
    QCOMPARE(guitest::replaceFolder("C:\\work\\case\\p.toml (save it)", "C:\\work\\case", "{case_dir}", true),
             std::string("{case_dir}/p.toml (save it)"));
    QCOMPARE(guitest::replaceFolder("/work/case/a\\b", "/work/case", "{case_dir}", false), std::string("{case_dir}/a\\b"));
#ifdef _WIN32
    const fs::path dir = core::pathFromUtf8("C:\\work\\case");
    QCOMPARE(guitest::normalise("out_dir = \"C:/work/case/MDATK_output\"\r\n", dir),
             std::string("out_dir = \"{case_dir}/MDATK_output\"\n"));
    QCOMPARE(guitest::normalise("Saved project C:\\work\\case\\sub\\p.mdatk.toml", dir),
             std::string("Saved project {case_dir}/sub/p.mdatk.toml"));
#else
    const fs::path dir = core::pathFromUtf8("/work/case");
    QCOMPARE(guitest::normalise("out_dir = \"/work/case/MDATK_output\"\r\n", dir),
             std::string("out_dir = \"{case_dir}/MDATK_output\"\n"));
#endif
}

void GuiTests::readsCaseFiles() {
    QTemporaryDir tmp;
    QVERIFY(tmp.isValid());
    const fs::path dir = core::pathFromUtf8(tmp.path().toStdString());
    const auto write = [&](const char* name, const char* text) {
        QFile f(QString::fromStdString(core::pathToUtf8(dir / name)));
        QVERIFY(f.open(QIODevice::WriteOnly));
        f.write(text);
    };
    write("bad_format.json", R"({"format": "other", "version": 1, "actions": []})");
    write("bad_action.json", R"({"format": "mdbridge-gui-case", "version": 1, "actions": [{"do": "open"}]})");
    write("bad_clock.json", R"({"format": "mdbridge-gui-case", "version": 1, "clock": "2026-10-15 14:32", "actions": []})");
    write("good.json", R"({"format": "mdbridge-gui-case", "version": 1, "clock": "2026-12-31T23:59:59",
                           "files": {"p.mdatk.toml": "projects/default.mdatk.toml"},
                           "actions": [{"action": "open", "path": "p.mdatk.toml"}, {"action": "save"}],
                           "expected": {"p.mdatk.toml": "gui/x.mdatk.toml"}})");
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::loadCase(dir / "bad_format.json", dir); })),
             QStringLiteral("bad_format: 'format' must be \"mdbridge-gui-case\""));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::loadCase(dir / "bad_action.json", dir); })),
             QStringLiteral("bad_action: action 1 must be an object with an \"action\" string"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::loadCase(dir / "bad_clock.json", dir); })),
             QStringLiteral("bad_clock: 'clock' must be YYYY-MM-DDTHH:MM:SS, not '2026-10-15 14:32'"));
    const guitest::GuiCase c = guitest::loadCase(dir / "good.json", dir);
    QCOMPARE(c.name, std::string("good"));
    QCOMPARE(c.actions.size(), std::size_t{2});
    QCOMPARE(c.files.size(), std::size_t{1});
    QCOMPARE(c.expected.front().second, std::string("gui/x.mdatk.toml"));
    QVERIFY(c.clock.has_value());
    {
        guitest::FixedClock clock(*c.clock);
        QCOMPARE(core::clockNow().isoformat(), std::string("2026-12-31T23:59:59"));
        QCOMPARE(core::clockNow().isoformat(), std::string("2027-01-01T00:00:00"));
    }
}

void GuiTests::cases_data() {
    QTest::addColumn<QString>("file");
    const fs::path dir = kData / "gui";
    const std::vector<fs::path> files = guitest::findCases(dir);
    QVERIFY2(!files.empty(), "no GUI cases in tests/data/gui");
    for (const fs::path& file : files) {
        const std::u8string name = file.lexically_relative(dir).replace_extension().generic_u8string();
        QTest::newRow(std::string(name.begin(), name.end()).c_str()) << QString::fromStdString(core::pathToUtf8(file));
    }
}

// The environment of a case (tests/data/README.md "gui"): XDG_CONFIG_HOME (and APPDATA) a folder of
// its own, HOME an empty folder, the case folder current; restored afterwards.
class CaseEnvironment {
public:
    CaseEnvironment(const fs::path& caseDir, const fs::path& configDir, const fs::path& homeDir) {
        for (const char* name : {"XDG_CONFIG_HOME", "APPDATA", "HOME", "USERPROFILE"}) {
            m_saved.emplace_back(name, qEnvironmentVariableIsSet(name) ? std::optional<QByteArray>(qgetenv(name)) : std::nullopt);
        }
        const QByteArray config = QByteArray::fromStdString(core::pathToUtf8(configDir));
        const QByteArray home = QByteArray::fromStdString(core::pathToUtf8(homeDir));
        qputenv("XDG_CONFIG_HOME", config);
        qputenv("APPDATA", config);
        qputenv("HOME", home);
        qputenv("USERPROFILE", home);
        m_cwd = fs::current_path();
        fs::current_path(caseDir);
    }
    ~CaseEnvironment() {
        std::error_code ec;
        fs::current_path(m_cwd, ec);
        for (const auto& [name, value] : m_saved) {
            if (value) qputenv(name, *value);
            else qunsetenv(name);
        }
    }
    CaseEnvironment(const CaseEnvironment&) = delete;
    CaseEnvironment& operator=(const CaseEnvironment&) = delete;

private:
    std::vector<std::pair<const char*, std::optional<QByteArray>>> m_saved;
    fs::path m_cwd;
};

std::string replaced(std::string text, const std::string& from, const std::string& to) {
    for (std::size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
    return text;
}

std::string readText(const fs::path& path) {
    QFile f(QString::fromStdString(core::pathToUtf8(path)));
    if (!f.open(QIODevice::ReadOnly)) throw CaseError("cannot read " + core::pathToUtf8(path));
    return f.readAll().toStdString();
}

void writeText(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    QFile f(QString::fromStdString(core::pathToUtf8(path)));
    if (!f.open(QIODevice::WriteOnly)) throw CaseError("cannot write " + core::pathToUtf8(path));
    f.write(text.data(), static_cast<qint64>(text.size()));
}

void GuiTests::cases() {
    QFETCH(QString, file);
    const guitest::GuiCase c = guitest::loadCase(core::pathFromUtf8(file.toStdString()), kData / "gui");
    if (const std::vector<std::string> missing = guitest::unsupportedNeeds(c); !missing.empty()) {
        std::string why = "needs what the harness does not have yet:";
        for (const std::string& need : missing) {
            why += " " + need;
        }
        QSKIP(why.c_str());
    }
    QTemporaryDir tmp;
    QVERIFY(tmp.isValid());
    const fs::path root = core::pathFromUtf8(tmp.path().toStdString());
    const fs::path caseDir = root / "case";
    const fs::path configDir = root / "config";
    const fs::path homeDir = root / "home";
    for (const fs::path& dir : {caseDir, configDir, homeDir}) fs::create_directories(dir);
    std::optional<guitest::FixedClock> clock;
    if (c.clock) clock.emplace(*c.clock);

    std::string failure;
    {
        CaseEnvironment environment(caseDir, configDir, homeDir);
        std::vector<std::string> log;
        std::optional<guitest::CaseRunner> runner;
        try {
            guitest::prepareCase(c, kData, caseDir);
            if (!c.settings.empty()) {
                // the generic form: a Windows path's '\\' would be an escape in the TOML string
                const std::u8string folder = caseDir.lexically_normal().generic_u8string();
                writeText(core::userSettingsPath(),
                          replaced(readText(kData / core::pathFromUtf8(c.settings)), "{case_dir}",
                                   std::string(folder.begin(), folder.end())));
            }
        } catch (const std::exception& e) {
            failure = e.what();
        }
        if (failure.empty()) {
            app::MainWindow window(render::Renderer::Software);
            window.resize(1320, 920);
            if (c.canvas) {  // the figure canvas at MDATK's size (its geometry, the map drags' pixels -> Å)
                window.preview()->setFixedSize(c.canvas->first, c.canvas->second);
                window.resize(std::max(1320, c.canvas->first + 600), std::max(920, c.canvas->second + 420));
            }
            window.show();
            QVERIFY(QTest::qWaitForWindowExposed(&window));
            try {
                runner.emplace(window, caseDir, screenshots(), c.name);
                runner->setOpenMayFail(c.needsAny("open_failure"));
                if (c.outcome.isTable()) {
                    if (const core::Value* records = c.outcome.find("canvas")) runner->setCanvasRecords(*records);
                }
                if (!c.start.empty()) runner->start(c.start);
                runner->run(c.actions);
                for (const std::string& note : runner->geometryNotes()) qInfo("%s", note.c_str());
                for (const fs::path& shot : runner->screenshots()) {
                    qInfo("screenshot %s", core::pathToUtf8(shot).c_str());
                }
            } catch (const std::exception& e) {
                failure = e.what();
            }
            window.logView()->flush();
            const QStringList lines = window.logView()->text().split(QLatin1Char('\n'), Qt::SkipEmptyParts);
            for (const QString& line : lines) log.push_back(line.toStdString());
            window.cancelLoad();
            if (failure.empty()) {
                failure = guitest::checkOutcome(c, *runner, log, kData, caseDir, configDir);
            }
            runner.reset();
        }
    }
    if (failure.empty()) {
        failure = guitest::checkExpected(c, kData, caseDir);
    }
    if (!failure.empty()) {
        tmp.setAutoRemove(false);
        failure += "\n(case folder kept: " + core::pathToUtf8(caseDir) + ")";
    }
    QVERIFY2(failure.empty(), failure.c_str());
}

int main(int argc, char** argv) {
    // offscreen: no display on Linux CI, no desktop on Windows CI; MDBRIDGE_GUI_TEST_PLATFORM
    // chooses another (e.g. xcb to watch the window)
    const QByteArray platform = qgetenv("MDBRIDGE_GUI_TEST_PLATFORM");
    qputenv("QT_QPA_PLATFORM", platform.isEmpty() ? QByteArray("offscreen") : platform);
    QLocale::setDefault(QLocale::c());  // numbers typed as the project file writes them
    QApplication app(argc, argv);
    guitest::useBundledFont();
    GuiTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_gui.moc"
