// mdbridge_gui_tests: the main window driven by action scripts (tests/data/gui, format in
// tests/data/README.md "gui"; tests/parity/make_gui_harness_cases.py and the per-page generators),
// the saved files compared with MDATK's, screenshots for CI artifacts. Qt Test on the "offscreen"
// platform (Linux and Windows, no display), the software 3D view.

#include <cstdlib>
#include <exception>
#include <memory>
#include <optional>
#include <string>

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

#include "gui/main_window.hpp"
#include "gui_harness.hpp"
#include "mdbridge/core/paths.hpp"

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
    combo->addItems({QStringLiteral("ps"), QStringLiteral("ns"), QStringLiteral("us")});
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
    QCOMPARE(combo->currentText(), QStringLiteral("ns"));
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
             QStringLiteral("inputs.time_unit: the combo box has no item 'fs'"));
    QCOMPARE(QString::fromStdString(errorOf([&] { guitest::setField(check, 1); })),
             QStringLiteral("video.enabled: a check box takes true or false"));
}

void GuiTests::comparesTexts() {
    QCOMPARE(guitest::compareText("a\nb\n", "a\nb\n"), std::string());
    QCOMPARE(guitest::compareText("a\nc\n", "a\nb\n"), std::string("line 2: expected 'b', got 'c'"));
    QCOMPARE(guitest::compareText("a\n", "a\nb\n"), std::string("line 2: expected 'b', got the end of the file"));
    QCOMPARE(guitest::compareText("a", "a\n"), std::string("the texts differ in their final line end"));
#ifdef _WIN32
    const fs::path dir = core::pathFromUtf8("C:\\work\\case");
    QCOMPARE(guitest::normalise("out_dir = \"C:/work/case/MDATK_output\"\r\n", dir),
             std::string("out_dir = \"{case_dir}/MDATK_output\"\n"));
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

void GuiTests::cases() {
    QFETCH(QString, file);
    const guitest::GuiCase c = guitest::loadCase(core::pathFromUtf8(file.toStdString()), kData / "gui");
    QTemporaryDir tmp;
    QVERIFY(tmp.isValid());
    const fs::path caseDir = core::pathFromUtf8(tmp.path().toStdString()) / "case";
    fs::create_directories(caseDir);
    std::optional<guitest::FixedClock> clock;
    if (c.clock) clock.emplace(*c.clock);

    std::string failure;
    {
        app::MainWindow window(render::Renderer::Software);
        window.resize(1280, 800);
        window.show();
        QVERIFY(QTest::qWaitForWindowExposed(&window));
        try {
            guitest::prepareCase(c, kData, caseDir);
            guitest::CaseRunner runner(window, caseDir, screenshots(), c.name);
            runner.run(c.actions);
            for (const fs::path& shot : runner.screenshots()) {
                qInfo("screenshot %s", core::pathToUtf8(shot).c_str());
            }
        } catch (const std::exception& e) {
            failure = e.what();
        }
        window.cancelLoad();
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
    GuiTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_gui.moc"
