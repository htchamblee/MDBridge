#pragma once
// The questions and messages of the main window (milestone 5 W1-A), behind one object the GUI tests
// replace: no modal dialog can be answered on the offscreen platform, and the action scripts
// (tests/data/README.md "gui") name the answers ("overwrite": "replace", "unsaved": "cancel" ...).
// File dialogs go through forms::runFileDialog (field_widgets.hpp), which tests replace the same way.
//
// The default Prompter shows Qt message boxes: questions with one button per choice (the first
// choice is the default button, Esc picks "cancel" when there is one), messages with their detail
// selectable (MDATK 0.1.0's dialogtext; long details scroll instead of being cut at 3 000
// characters as in v0.0.4.1, AP-18), and prints to stderr.

#include <memory>
#include <utility>
#include <vector>

#include <QString>

class QMessageBox;
class QWidget;

namespace mdbridge::app {

struct Question {
    // What is asked (the GUI cases' "answers" keys): "unsaved" (Save / Discard / Cancel before New,
    // Open, Close project and closing the window: AP-13), "overwrite" (Save as onto a file that
    // exists: AP-11), "external" (Save onto a file another program changed), "reload", "recover"
    // (a recovery copy found at the start: D6).
    QString kind;
    QString path;   // the file the question is about ("" for a project without one)
    QString title;  // the dialog's title
    QString text;   // the question
    QString detail;
    // (id, button label); the first is the answer when none is given (GUI tests)
    std::vector<std::pair<QString, QString>> choices;
    // The dialog's default button (Return); "" = the first choice. "cancel" for questions whose
    // first choice replaces or discards something.
    QString defaultId;
};

enum class MessageKind { Error, Info };

class Prompter {
public:
    Prompter() = default;
    virtual ~Prompter() = default;
    Prompter(const Prompter&) = delete;
    Prompter& operator=(const Prompter&) = delete;

    // The id of the choice taken.
    virtual QString ask(QWidget* parent, const Question& question);
    // A modal message (the window has already logged it).
    virtual void message(QWidget* parent, MessageKind kind, const QString& text, const QString& detail);
    // A line for the terminal (MDATK's "warning: could not load ..." at the start).
    virtual void printError(const QString& line);
};

// The Prompter windows use unless they are given another (the Qt dialogs above).
Prompter& defaultPrompter();

// The message box Prompter::message shows, not yet shown (the tests' screenshots of the result
// dialogs show it with show() instead of exec()).
std::unique_ptr<QMessageBox> messageBox(QWidget* parent, MessageKind kind, const QString& text, const QString& detail);

}  // namespace mdbridge::app
