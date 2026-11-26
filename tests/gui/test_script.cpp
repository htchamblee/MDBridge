// mdbridge_gui_script_tests: `mdbridge gui --script` (src/app/gui/script_runner.*), the action
// scripts the smoke test of finished builds runs (packaging/smoke_test.py --gui): a value typed into
// a field and saved changes exactly that line of the file, every page is photographed, a failing
// action stops the script with exit code 1 and names the action, a malformed script is refused.
// Qt Test on the offscreen platform.

#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <QApplication>
#include <QLocale>
#include <QTemporaryDir>
#include <QTest>

#include "gui/display_path.hpp"
#include "gui/main_window.hpp"
#include "gui/script_runner.hpp"
#include "mdbridge/core/paths.hpp"

namespace fs = std::filesystem;
using namespace mdbridge;

namespace {

void writeFile(const fs::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out << text;
}

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

std::vector<std::string> lines(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();  // project files are CRLF on Windows
        out.push_back(line);
    }
    return out;
}

}  // namespace

class ScriptTests : public QObject {
    Q_OBJECT

private slots:
    void initTestCase();
    void editSaveAndScreenshots();
    void failingActionStops();
    void malformedScriptRefused();
    void problemShowsPathsAsTheInterface();

private:
    QTemporaryDir m_config;
};

void ScriptTests::initTestCase() {
    QVERIFY(m_config.isValid());
    qputenv("XDG_CONFIG_HOME", m_config.path().toUtf8());
    qputenv("APPDATA", m_config.path().toUtf8());
}

void ScriptTests::editSaveAndScreenshots() {
    QTemporaryDir dir;
    const fs::path folder = core::pathFromUtf8(dir.path().toStdString());
    const fs::path project = folder / "p.mdatk.toml";
    writeFile(project, "format = \"mdatk-project\"\n");
    writeFile(folder / "s.json", R"({"format": "mdbridge-gui-script", "version": 1, "actions": [
        {"action": "wait"},
        {"action": "save"},
        {"action": "copy_project", "to": "before.toml"},
        {"action": "set", "key": "output.video.fps", "value": 30},
        {"action": "save"},
        {"action": "screenshot", "name": "page", "pages": true}]})");
    app::MainWindow window(render::Renderer::Software);
    QVERIFY(window.startWith(project));
    window.show();
    app::ScriptRunner runner(window, folder / "s.json", folder / "shots");
    QCOMPARE(runner.run(), 0);

    const std::vector<std::string> before = lines(readFile(folder / "before.toml"));
    const std::vector<std::string> after = lines(readFile(project));
    QCOMPARE(before.size(), after.size());
    std::vector<std::string> changed;
    for (std::size_t k = 0; k < before.size(); ++k) {
        if (before[k] != after[k] && before[k].rfind("saved = ", 0) != 0) changed.push_back(before[k] + " -> " + after[k]);
    }
    QCOMPARE(changed, std::vector<std::string>{"fps = 24 -> fps = 30"});
    QVERIFY(!window.document().isModified());
    for (const app::MainWindow::Page& page : window.pages()) {
        const fs::path shot = folder / "shots" / ("page-" + page.id.toStdString() + ".png");
        QVERIFY2(fs::is_regular_file(shot) && fs::file_size(shot) > 1000, shot.string().c_str());
    }
}

void ScriptTests::failingActionStops() {
    QTemporaryDir dir;
    const fs::path folder = core::pathFromUtf8(dir.path().toStdString());
    const fs::path project = folder / "p.mdatk.toml";
    writeFile(project, "format = \"mdatk-project\"\n");
    writeFile(folder / "s.json", R"({"actions": [
        {"action": "set", "key": "no.such.field", "value": 1},
        {"action": "save"}]})");
    const std::string original = readFile(project);
    app::MainWindow window(render::Renderer::Software);
    QVERIFY(window.startWith(project));
    window.show();
    app::ScriptRunner runner(window, folder / "s.json", {});
    QCOMPARE(runner.run(), 1);
    QCOMPARE(readFile(project), original);  // the save after the failing action did not run
}

void ScriptTests::malformedScriptRefused() {
    QTemporaryDir dir;
    const fs::path folder = core::pathFromUtf8(dir.path().toStdString());
    writeFile(folder / "s.json", R"({"steps": []})");
    app::MainWindow window(render::Renderer::Software);
    QVERIFY_THROWS_EXCEPTION(std::runtime_error, app::ScriptRunner(window, folder / "s.json", {}));
    writeFile(folder / "t.json", "{");
    QVERIFY_THROWS_EXCEPTION(std::exception, app::ScriptRunner(window, folder / "t.json", {}));
}

// The report's "problem" (and stderr) shows a path as the interface does: gui::displayPath, '\'
// throughout on Windows (WN-12), unchanged on Linux and macOS.
void ScriptTests::problemShowsPathsAsTheInterface() {
    QTemporaryDir dir;
    const fs::path folder = core::pathFromUtf8(dir.path().toStdString());
    const fs::path project = folder / "p.mdatk.toml";
    writeFile(project, "format = \"mdatk-project\"\n");
    writeFile(folder / "s.json", R"({"actions": [{"action": "copy_project", "to": "missing/copy.toml"}]})");
    app::MainWindow window(render::Renderer::Software);
    QVERIFY(window.startWith(project));
    app::ScriptRunner runner(window, folder / "s.json", {});
    QCOMPARE(runner.run(), 1);
    const fs::path target = fs::absolute(folder / "s.json").parent_path() / core::pathFromUtf8("missing/copy.toml");
    const std::string expected = app::gui::displayPath(target).toStdString();
    const std::string& problem = runner.problem();
    QVERIFY2(problem.find("could not copy the project file to " + expected + ": ") != std::string::npos, problem.c_str());
#ifdef _WIN32
    QVERIFY2(problem.find('/') == std::string::npos, problem.c_str());
#endif
}

int main(int argc, char** argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QLocale::setDefault(QLocale::c());
    QApplication app(argc, argv);
    QApplication::setStyle(QStringLiteral("Fusion"));
    ScriptTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "test_script.moc"
