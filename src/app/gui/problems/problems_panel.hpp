#pragma once
// Where the window shows problems (milestone 5 W1-D; MDATK 0.1.1's gui/problemview.py, the 0.0.10
// presentation of owner decision D7): a count badge per page button, and the problems list, which
// names each problem's place ("Page › Section › Field") and message and jumps to the field when one
// is clicked. What the problems are is FormContext::allProblems() (the checks of validate(), with
// v0.0.4.1's blocking rules, checked as you type, and the field widgets' own warnings); the field
// widgets flag themselves.
//
// The jump goes through a reveal function the window gives (ProblemsPanel::setRevealFunction): it
// shows the page of the field (core::problemPage) and calls that page's revealField(field), which
// can use revealFieldIn() below.

#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include <QLabel>
#include <QString>
#include <QTimer>
#include <QWidget>

#include "gui/forms/form_context.hpp"
#include "mdbridge/core/problems.hpp"

class QListWidget;
class QListWidgetItem;

namespace mdbridge::app::problems {

// field_parents of MDATK 0.1.1: "groups.1.residues" -> "groups.1.residues", "groups.1", "groups".
std::vector<std::string> fieldParents(std::string_view field);

struct ProblemCount {
    int errors = 0;
    int warnings = 0;
    int total() const { return errors + warnings; }
    // "error" when any is an error, else "warning" (MDATK 0.1.1's _level).
    QString level() const { return errors > 0 ? QStringLiteral("error") : QStringLiteral("warning"); }
};

// The problems a page shows (Problem::page() == page).
ProblemCount countOn(const std::vector<core::Problem>& problems, std::string_view page);
std::vector<core::Problem> problemsOn(const std::vector<core::Problem>& problems, std::string_view page);
// "2 errors, 1 warning" ("" for none).
QString countText(const ProblemCount& count);
// The lines a page button's tooltip adds (MDATK 0.1.1's _summary): "• <location without the page>:
// <message>" for the first `limit`, then "… and <n> more".
QString summaryText(const std::vector<core::Problem>& problems, int limit = 6);

// Shows the widget of `field` below `root` (a list cell or field widget named after the path; else
// that of the nearest parent path), scrolling every scroll area and stacked widget above it so that
// it is in view, and gives it the keyboard focus (`focus`; a line edit's text is not selected).
// False when nothing below `root` edits the field or a parent of it. Pages use it for revealField.
bool revealFieldIn(QWidget* root, const QString& field, bool focus = true);

// A page button's count badge (MDATK 0.1.1's badge_overlay / set_badge): the number of the page's
// problems ("99+" from 100), red when one is an error, amber for warnings only, hidden at 0. Its
// tooltip is summaryText() of the page's problems. objectName "badge:<page>".
class ProblemBadge : public QLabel {
    Q_OBJECT

public:
    ProblemBadge(forms::FormContext& context, QString page, QWidget* parent = nullptr);
    const QString& page() const { return m_page; }
    const ProblemCount& count() const { return m_count; }
    // Counts again now (it does so once the event loop runs after the problems change).
    void refresh();

signals:
    void countChanged(int errors, int warnings);

private:
    forms::FormContext& m_context;
    QString m_page;
    ProblemCount m_count;
    QTimer m_timer;
};

// The problems list (MDATK 0.1.1's problems button and popover, as a panel the window places): a
// summary line ("2 errors, 1 warning" / "No problems"), then one row per problem in
// FormContext::allProblems()'s order with the level's icon, the location in bold ("Project" for a
// problem no setting names) and the message; clicking a row (or Enter on it) reveals its field.
// objectNames: "problems", "problems:summary", "problems:list".
class ProblemsPanel : public QWidget {
    Q_OBJECT

public:
    // Roles of a row's data: the field path, the page id and the level.
    static constexpr int FieldRole = Qt::UserRole + 1;
    static constexpr int PageRole = Qt::UserRole + 2;
    static constexpr int LevelRole = Qt::UserRole + 3;

    explicit ProblemsPanel(forms::FormContext& context, QWidget* parent = nullptr);

    // How a click shows a field: the window shows the field's page and calls its revealField.
    using RevealFunction = std::function<bool(const QString& field)>;
    void setRevealFunction(RevealFunction reveal) { m_reveal = std::move(reveal); }
    // Shows the field of a problem (the reveal function; false without one or when it finds none).
    bool revealProblem(const QString& field);

    // The problems shown, in the list's order.
    const std::vector<core::Problem>& problems() const { return m_problems; }
    ProblemCount count() const;
    ProblemCount pageCount(const QString& page) const;
    // The tooltip lines of a page button (summaryText of the page's problems; "" when none).
    QString pageSummary(const QString& page) const;
    // The next problem of a page a click on its badge shows (MDATK 0.1.1 cycles through them);
    // nullptr when the page has none.
    const core::Problem* nextOnPage(const QString& page);

    QListWidget* list() const { return m_list; }
    QLabel* summaryLabel() const { return m_summary; }
    // Reads the problems again now (it does so once the event loop runs after they change).
    void refresh();

signals:
    // After each refresh: the counts (the window's page buttons and status may show them).
    void problemsShown(int errors, int warnings);
    // A row was activated (before the reveal function runs).
    void problemActivated(const QString& field, const QString& page);

private:
    void activate(QListWidgetItem* item);

    forms::FormContext& m_context;
    RevealFunction m_reveal;
    std::vector<core::Problem> m_problems;
    std::map<QString, int> m_next;
    QLabel* m_summary = nullptr;
    QListWidget* m_list = nullptr;
    QTimer m_timer;
};

}  // namespace mdbridge::app::problems
