#include "script_runner.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <QAbstractButton>
#include <QAbstractSpinBox>
#include <QApplication>
#include <QComboBox>
#include <QDeadlineTimer>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPixmap>
#include <QPlainTextEdit>
#include <QThread>
#include <QWidget>

#include "gui/actions/window_actions.hpp"
#include "gui/dialogs/image_dialog.hpp"
#include "gui/dialogs/prompts.hpp"
#include "gui/display_path.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "gui/jobs/job_runner.hpp"
#include "gui/main_window.hpp"
#include "gui/preview_widget.hpp"
#include "loaded_libraries.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"

namespace mdbridge::app {

namespace fs = std::filesystem;

namespace {

constexpr int kLoadTimeoutMs = 180000;
constexpr int kStillTimeoutMs = 600000;

// An action failed: the message says why.
class ActionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

QString qs(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }
std::string utf8(const QString& text) { return forms::toStdString(text); }
// A path in a message (the report's "problem", stderr) as the interface shows it: native
// separators on Windows (WN-12); the report's path fields keep the stored form.
std::string displayed(const fs::path& path) { return utf8(gui::displayPath(path)); }

// Runs the event loop until `done` holds or `timeoutMs` passed (as QTest::qWaitFor).
bool waitFor(const std::function<bool()>& done, int timeoutMs) {
    QDeadlineTimer deadline(timeoutMs);
    while (!done()) {
        if (deadline.hasExpired()) return false;
        QCoreApplication::processEvents(QEventLoop::AllEvents, 20);
        QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
        QThread::msleep(5);
    }
    return true;
}

void settle() {
    QApplication::processEvents();
    QApplication::processEvents();
}

// A key pressed and released on `widget`, as the keyboard sends it.
void key(QWidget* widget, int code, const QString& text = QString()) {
    QKeyEvent press(QEvent::KeyPress, code, Qt::NoModifier, text);
    QApplication::sendEvent(widget, &press);
    QKeyEvent release(QEvent::KeyRelease, code, Qt::NoModifier, text);
    QApplication::sendEvent(widget, &release);
}

// `text` typed into `widget`, one key per character.
void typeKeys(QWidget* widget, const QString& text) {
    for (const QChar c : text) {
        const char16_t u = c.toUpper().unicode();
        key(widget, u >= 0x20 && u < 0x7f ? static_cast<int>(u) : static_cast<int>(Qt::Key_unknown), QString(c));
    }
}

// What a user types for `value`: strings as they are, numbers as the project file writes them,
// lists of strings comma-separated.
QString typedText(const core::Value& value) {
    if (const auto* s = value.get<std::string>()) return qs(*s);
    if (const auto* i = value.get<std::int64_t>()) return QString::number(*i);
    if (const auto* d = value.get<double>()) return qs(core::pythonFloatRepr(*d));
    if (const auto* items = value.get<core::Array>()) {
        QStringList parts;
        for (const core::Value& item : *items) parts.append(typedText(item));
        return parts.join(QStringLiteral(", "));
    }
    throw ActionError("cannot type a " + std::string(core::typeName(value.type())));
}

// Select all, type, Enter: line edits, spin boxes, editable drop-downs.
void typeLine(QWidget* widget, const QString& text) {
    QWidget* target = widget;
    if (auto* box = qobject_cast<QComboBox*>(widget)) {
        target = box->lineEdit();
        box->lineEdit()->selectAll();
    } else if (auto* spin = qobject_cast<QAbstractSpinBox*>(widget)) {
        spin->selectAll();
    } else if (auto* edit = qobject_cast<QLineEdit*>(widget)) {
        edit->selectAll();
    }
    if (text.isEmpty()) {
        key(target, Qt::Key_Delete);
    } else {
        typeKeys(target, text);
    }
    key(target, Qt::Key_Return, QStringLiteral("\r"));
}

// Multi-line text is stored when the field is left: as its field widget does then.
void commitAsOnLeaving(QWidget* widget) {
    for (QWidget* w = widget; w; w = w->parentWidget()) {
        if (auto* field = qobject_cast<forms::FieldWidget*>(w)) {
            field->commit();
            return;
        }
    }
    widget->clearFocus();
}

void enterValue(QWidget* widget, const core::Value& value) {
    const std::string name = utf8(widget->objectName());
    widget->setFocus(Qt::OtherFocusReason);
    if (auto* box = qobject_cast<QComboBox*>(widget)) {
        const QString text = typedText(value);
        if (box->isEditable()) {
            typeLine(box, text);
            return;
        }
        int index = box->findData(text);
        if (index < 0) index = box->findText(text);
        if (index < 0) throw ActionError(name + ": the drop-down has no item '" + utf8(text) + "'");
        box->setCurrentIndex(index);
        emit box->activated(index);  // as when a user picks the item
        return;
    }
    if (qobject_cast<QLineEdit*>(widget) || qobject_cast<QAbstractSpinBox*>(widget)) {
        typeLine(widget, typedText(value));
        return;
    }
    if (auto* plain = qobject_cast<QPlainTextEdit*>(widget)) {
        plain->selectAll();
        const QStringList lines = typedText(value).split(QLatin1Char('\n'));
        for (qsizetype k = 0; k < lines.size(); ++k) {
            if (k > 0) key(plain, Qt::Key_Return, QStringLiteral("\r"));
            typeKeys(plain, lines[k]);
        }
        commitAsOnLeaving(plain);
        return;
    }
    if (auto* button = qobject_cast<QAbstractButton*>(widget); button && button->isCheckable()) {
        const bool* wanted = value.get<bool>();
        if (!wanted) throw ActionError(name + ": a check box takes true or false");
        if (button->isChecked() != *wanted) key(button, Qt::Key_Space, QStringLiteral(" "));
        if (button->isChecked() != *wanted) throw ActionError(name + ": Space did not change the check box");
        return;
    }
    throw ActionError(name + ": cannot enter a value into a " + widget->metaObject()->className());
}

const std::string& stringEntry(const core::Value& action, std::string_view k) {
    const core::Value* v = action.find(k);
    const auto* s = v ? v->get<std::string>() : nullptr;
    if (!s) throw ActionError("'" + std::string(k) + "' must be a string");
    return *s;
}

int intEntry(const core::Value& action, std::string_view k, int fallback) {
    const core::Value* v = action.find(k);
    if (!v) return fallback;
    if (const auto* i = v->get<std::int64_t>()) return static_cast<int>(*i);
    if (const auto* d = v->get<double>()) return static_cast<int>(*d);
    throw ActionError("'" + std::string(k) + "' must be a number");
}

std::string readText(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + displayed(path));
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// The questions and messages of the window while the script runs: answered and recorded, never a
// modal dialog (nobody is there to answer it).
class ScriptPrompter : public Prompter {
public:
    QString ask(QWidget*, const Question& question) override {
        QString answer;
        if (const core::Value* a = answers.isNull() ? nullptr : answers.find(utf8(question.kind))) {
            answer = qs(a->stringOr());
        }
        if (answer.isEmpty()) answer = question.defaultId;
        if (answer.isEmpty() && !question.choices.empty()) answer = question.choices.front().first;
        questions.push_back(core::Table{{"kind", utf8(question.kind)},
                                        {"text", utf8(question.text)},
                                        {"answer", utf8(answer)}});
        return answer;
    }
    void message(QWidget*, MessageKind kind, const QString& text, const QString& detail) override {
        messages.push_back(core::Table{{"kind", kind == MessageKind::Error ? "error" : "info"},
                                       {"text", utf8(text)},
                                       {"detail", utf8(detail)}});
        if (kind == MessageKind::Error) ++errors;
    }
    void printError(const QString& line) override { std::cerr << utf8(line) << '\n'; }

    core::Value answers;
    core::Array questions;
    core::Array messages;
    int errors = 0;
};

}  // namespace

struct ScriptRunner::Impl {
    MainWindow& window;
    fs::path scriptDir;
    fs::path output;
    core::Value script;
    ScriptPrompter prompter;
    core::Array results;
    std::string lastProblem;  // the report's "problem" of the last run

    Impl(MainWindow& w, const fs::path& file, const fs::path& out) : window(w) {
        script = core::parseJson(readText(file), core::pathToUtf8(file));
        scriptDir = file.parent_path().empty() ? fs::current_path() : fs::absolute(file).parent_path();
        output = out.empty() ? scriptDir : fs::absolute(out);
        const core::Value* actions = script.find("actions");
        if (!actions || !actions->get<core::Array>()) {
            throw std::runtime_error(displayed(file) + ": the script has no \"actions\" list");
        }
    }

    fs::path resolved(const std::string& path) const {
        const fs::path p = core::pathFromUtf8(path);
        return p.is_absolute() ? p : scriptDir / p;
    }

    void waitUntilLoaded(int timeoutMs) {
        if (!waitFor([this] { return !window.loading() && window.previewSettled(); }, timeoutMs)) {
            throw ActionError("the project and its preview were not ready in " + std::to_string(timeoutMs / 1000) + " s");
        }
    }

    fs::path grab(QWidget* widget, const std::string& name) {
        settle();
        std::error_code ec;
        fs::create_directories(output, ec);
        const fs::path file = output / core::pathFromUtf8(name + ".png");
        const QPixmap pixmap = widget->grab();
        if (pixmap.isNull() || !pixmap.save(qs(core::pathToUtf8(file)), "PNG")) {
            throw ActionError("could not save the screenshot " + displayed(file));
        }
        return file;
    }

    // Runs one action; `result` gets what it did.
    void runAction(const core::Value& action, core::Table& result) {
        const std::string& what = stringEntry(action, "action");
        if (what == "wait") {
            waitUntilLoaded(intEntry(action, "timeout_s", kLoadTimeoutMs / 1000) * 1000);
        } else if (what == "page") {
            const std::string& page = stringEntry(action, "page");
            if (!window.showPage(qs(page))) throw ActionError("there is no page '" + page + "'");
        } else if (what == "set" || what == "type") {
            const QString k = qs(stringEntry(action, "key"));
            window.revealField(k);
            settle();
            QList<QWidget*> found = window.findChildren<QWidget*>(k);
            if (found.size() != 1) {
                throw ActionError(std::to_string(found.size()) + " widgets are named '" + utf8(k) + "'");
            }
            if (what == "set") {
                const core::Value* value = action.find("value");
                if (!value) throw ActionError("'value' is missing");
                enterValue(found.front(), *value);
            } else {
                found.front()->setFocus(Qt::OtherFocusReason);
                typeLine(found.front(), qs(stringEntry(action, "text")));
            }
            settle();
            const core::Value stored = window.document().get(utf8(k));
            result.emplace_back("stored", stored);
        } else if (what == "save") {
            if (!window.save()) throw ActionError("Save did not save the project");
            result.emplace_back("file", core::pathToUtf8(window.projectFile()));
        } else if (what == "copy_project") {
            const fs::path to = resolved(stringEntry(action, "to"));
            std::error_code ec;
            fs::copy_file(window.projectFile(), to, fs::copy_options::overwrite_existing, ec);
            if (ec) throw ActionError("could not copy the project file to " + displayed(to) + ": " + ec.message());
            result.emplace_back("file", core::pathToUtf8(to));
        } else if (what == "render_still") {
            renderStill(action, result);
        } else if (what == "screenshot") {
            const std::string name = action.find("name") ? stringEntry(action, "name") : std::string("window");
            core::Array files;
            const core::Value* pages = action.find("pages");
            if (pages && pages->get<bool>() && *pages->get<bool>()) {
                const QString before = window.currentPage();
                for (const MainWindow::Page& page : window.pages()) {
                    window.showPage(page.id);
                    waitUntilLoaded(kLoadTimeoutMs);
                    files.push_back(core::pathToUtf8(grab(&window, name + "-" + utf8(page.id))));
                }
                window.showPage(before);
            } else {
                files.push_back(core::pathToUtf8(grab(&window, name)));
            }
            result.emplace_back("files", std::move(files));
        } else {
            throw ActionError("unknown action (src/app/gui/script_runner.hpp lists the actions)");
        }
    }

    void renderStill(const core::Value& action, core::Table& result) {
        WindowActions* actions = window.windowActions();
        jobs::JobRunner* jobsRunner = window.jobRunner();
        QString shown;
        QString failure;
        const QMetaObject::Connection a =
            QObject::connect(actions, &WindowActions::stillShown, actions, [&shown](const QString& path) { shown = path; });
        const QMetaObject::Connection b = QObject::connect(
            jobsRunner, &jobs::JobRunner::errorReported, actions,
            [&failure](const QString& text, const QString& detail) { failure = text + QStringLiteral(": ") + detail; });
        const int errorsBefore = prompter.errors;
        const bool started = actions->renderStill(intEntry(action, "frame", 0));
        const bool done = started && waitFor([&] { return !shown.isEmpty() || !failure.isEmpty() || !jobsRunner->running(); },
                                             intEntry(action, "timeout_s", kStillTimeoutMs / 1000) * 1000);
        waitFor([&] { return !jobsRunner->running(); }, 10000);
        QObject::disconnect(a);
        QObject::disconnect(b);
        if (!started) throw ActionError("Render still was refused (see the messages)");
        if (!done) throw ActionError("the still was not rendered in time");
        if (!failure.isEmpty()) throw ActionError("Render still failed: " + utf8(failure));
        if (shown.isEmpty() || prompter.errors > errorsBefore) throw ActionError("Render still showed no image");
        result.emplace_back("image", utf8(shown));
        if (const core::Value* to = action.find("to")) {
            const fs::path target = resolved(to->stringOr());
            std::error_code ec;
            fs::copy_file(core::pathFromUtf8(utf8(shown)), target, fs::copy_options::overwrite_existing, ec);
            if (ec) throw ActionError("could not copy the still to " + displayed(target) + ": " + ec.message());
            result.emplace_back("file", core::pathToUtf8(target));
        }
        if (ImageDialog* dialog = actions->stillDialog()) {
            if (const core::Value* shot = action.find("screenshot")) {
                result.emplace_back("dialog", core::pathToUtf8(grab(dialog, shot->stringOr("still-dialog"))));
            }
            dialog->close();
        }
    }

    int run() {
        Prompter* previous = &window.prompter();
        window.setPrompter(&prompter);
        std::string problem;
        const core::Array& actions = *script.find("actions")->get<core::Array>();
        for (std::size_t k = 0; k < actions.size() && problem.empty(); ++k) {
            const core::Value& action = actions[k];
            const core::Value* answers = action.find("answers");
            prompter.answers = answers ? *answers : core::Value();
            const core::Value* name = action.find("action");
            core::Table result{{"action", name ? name->stringOr("?") : std::string("?")}};
            const int errorsBefore = prompter.errors;
            try {
                runAction(action, result);
                if (prompter.errors > errorsBefore) {
                    const core::Value& last = prompter.messages.back();
                    throw ActionError("the window reported an error: " + last.find("text")->stringOr() + " (" +
                                      last.find("detail")->stringOr() + ")");
                }
                result.emplace_back("ok", true);
            } catch (const std::exception& e) {
                // paths in the reasons other code gives (the pipeline's, the file system's) shown too
                const std::string reason = utf8(gui::displayPathsIn(qs(e.what())));
                problem = "action " + std::to_string(k + 1) + " (" + result.front().second.stringOr() + "): " + reason;
                result.emplace_back("ok", false);
                result.emplace_back("problem", reason);
            }
            results.push_back(std::move(result));
        }
        window.setPrompter(previous);
        core::Array libraries;
        for (const std::string& lib : loadedLibraries()) libraries.push_back(lib);
        core::Array pageIds;
        for (const MainWindow::Page& page : window.pages()) pageIds.push_back(utf8(page.id));
        const core::Value report = core::Table{
            {"ok", problem.empty()},
            {"problem", problem},
            {"project", core::pathToUtf8(window.projectFile())},
            {"modified", window.document().isModified()},
            {"pages", std::move(pageIds)},
            {"actions", results},
            {"messages", prompter.messages},
            {"questions", prompter.questions},
            {"libraries", std::move(libraries)},
        };
        lastProblem = problem;
        std::cout << core::toJson(report, -1) << std::endl;
        if (!problem.empty()) std::cerr << "mdbridge: script: " << problem << '\n';
        return problem.empty() ? 0 : 1;
    }
};

ScriptRunner::ScriptRunner(MainWindow& window, const fs::path& script, const fs::path& output)
    : m_impl(std::make_unique<Impl>(window, script, output)) {}

ScriptRunner::~ScriptRunner() = default;

int ScriptRunner::run() { return m_impl->run(); }

const std::string& ScriptRunner::problem() const { return m_impl->lastProblem; }

}  // namespace mdbridge::app
