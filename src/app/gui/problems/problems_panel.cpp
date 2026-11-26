// The problems list and page badges: MDATK 0.1.1's gui/problemview.py (owner's code, ported).

#include "gui/problems/problems_panel.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>

#include <QAbstractItemView>
#include <QFont>
#include <QItemSelectionModel>
#include <QListWidget>
#include <QListWidgetItem>
#include <QScrollArea>
#include <QStackedWidget>
#include <QStyle>
#include <QTableView>
#include <QVBoxLayout>

#include "gui/display_path.hpp"
#include "gui/forms/field_widgets.hpp"
#include "gui/forms/list_editor.hpp"

namespace mdbridge::app::problems {

using forms::toQString;
using forms::toStdString;

namespace {

// FLAG_COLORS["…"]["badge"] of MDATK 0.1.1
constexpr const char* kErrorBadge = "#e01b24";
constexpr const char* kWarningBadge = "#c64600";

std::optional<std::size_t> rowIndex(std::string_view text) {
    if (text.empty() || text.size() > 9) {
        return std::nullopt;
    }
    std::size_t n = 0;
    for (const char c : text) {
        if (c < '0' || c > '9') {
            return std::nullopt;
        }
        n = n * 10 + static_cast<std::size_t>(c - '0');
    }
    return n;
}

// Every scroll area and stacked widget above `widget` (up to `root`) shows it.
void bringIntoView(QWidget* root, QWidget* widget) {
    QWidget* child = widget;
    for (QWidget* w = widget->parentWidget(); w; w = w->parentWidget()) {
        if (auto* stack = qobject_cast<QStackedWidget*>(w)) {
            stack->setCurrentWidget(child);
        }
        child = w;
        if (w == root) {
            break;
        }
    }
    for (QWidget* w = widget->parentWidget(); w; w = w->parentWidget()) {
        if (auto* scroll = qobject_cast<QScrollArea*>(w)) {
            scroll->ensureWidgetVisible(widget, 20, 20);
        }
        if (w == root) {
            break;
        }
    }
}

QIcon levelIcon(const QWidget* widget, const std::string& level) {
    return widget->style()->standardIcon(level == "warning" ? QStyle::SP_MessageBoxWarning
                                                            : QStyle::SP_MessageBoxCritical);
}

}  // namespace

std::vector<std::string> fieldParents(std::string_view field) {
    std::vector<std::string> out;
    if (field.empty()) {
        return out;
    }
    std::string_view rest = field;
    while (true) {
        out.emplace_back(rest);
        const std::size_t dot = rest.rfind('.');
        if (dot == std::string_view::npos) {
            break;
        }
        rest = rest.substr(0, dot);
    }
    return out;
}

ProblemCount countOn(const std::vector<core::Problem>& problems, std::string_view page) {
    ProblemCount out;
    for (const core::Problem& p : problems) {
        if (p.page() == page) {
            (p.level == "warning" ? out.warnings : out.errors) += 1;
        }
    }
    return out;
}

std::vector<core::Problem> problemsOn(const std::vector<core::Problem>& problems, std::string_view page) {
    std::vector<core::Problem> out;
    for (const core::Problem& p : problems) {
        if (p.page() == page) {
            out.push_back(p);
        }
    }
    return out;
}

QString countText(const ProblemCount& count) {
    QStringList parts;
    if (count.errors) {
        parts << QStringLiteral("%1 error%2").arg(count.errors).arg(count.errors != 1 ? QStringLiteral("s") : QString());
    }
    if (count.warnings) {
        parts << QStringLiteral("%1 warning%2")
                     .arg(count.warnings)
                     .arg(count.warnings != 1 ? QStringLiteral("s") : QString());
    }
    return parts.join(QStringLiteral(", "));
}

QString summaryText(const std::vector<core::Problem>& problems, int limit) {
    QStringList lines;
    const std::size_t shown = std::min(problems.size(), static_cast<std::size_t>(std::max(limit, 0)));
    for (std::size_t i = 0; i < shown; ++i) {
        const QString location = toQString(problems[i].location);
        const qsizetype cut = location.indexOf(QStringLiteral(" › "));
        const QString where = cut < 0 ? location : location.mid(cut + 3);
        lines << QStringLiteral("• %1: %2").arg(where, gui::displayPathsIn(toQString(problems[i].message)));
    }
    if (problems.size() > shown) {
        lines << QStringLiteral("… and %1 more").arg(problems.size() - shown);
    }
    return lines.join(QLatin1Char('\n'));
}

bool revealFieldIn(QWidget* root, const QString& field, bool focus) {
    if (!root || field.isEmpty()) {
        return false;
    }
    const std::vector<forms::ListEditor*> lists = [&] {
        std::vector<forms::ListEditor*> out;
        for (forms::ListEditor* e : root->findChildren<forms::ListEditor*>()) {
            out.push_back(e);
        }
        return out;
    }();
    for (const std::string& path : fieldParents(toStdString(field))) {
        if (auto* widget = root->findChild<forms::FieldWidget*>(toQString(path) + QStringLiteral(":field"))) {
            bringIntoView(root, widget);
            if (focus) {
                widget->setFocus(Qt::OtherFocusReason);  // not Tab: a line edit's text stays unselected
            }
            return true;
        }
        const forms::CellLocation cell = forms::findCell(root, path);
        if (cell.editor) {
            bringIntoView(root, cell.editor);
            QTableView* view = cell.editor->view();
            view->setCurrentIndex(cell.index);
            view->scrollTo(cell.index);
            if (focus) {
                view->setFocus(Qt::OtherFocusReason);
            }
            return true;
        }
        // a row of a list editor ("groups.1")
        for (forms::ListEditor* editor : lists) {
            const std::string& list = editor->path();
            if (path.size() > list.size() + 1 && path.compare(0, list.size(), list) == 0 && path[list.size()] == '.') {
                const auto row = rowIndex(std::string_view(path).substr(list.size() + 1));
                if (row && *row < static_cast<std::size_t>(editor->rowCount())) {
                    bringIntoView(root, editor);
                    editor->selectRow(static_cast<int>(*row));
                    if (focus) {
                        editor->view()->setFocus(Qt::OtherFocusReason);
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

// ---- ProblemBadge ----------------------------------------------------------------------------

ProblemBadge::ProblemBadge(forms::FormContext& context, QString page, QWidget* parent)
    : QLabel(parent), m_context(context), m_page(std::move(page)) {
    setObjectName(QStringLiteral("badge:") + m_page);
    setAlignment(Qt::AlignCenter);
    setAttribute(Qt::WA_TransparentForMouseEvents);  // the page button under it takes the click
    QFont f = font();
    f.setBold(true);
    f.setPointSizeF(std::max(6.0, f.pointSizeF() * 0.8));
    setFont(f);
    m_timer.setSingleShot(true);
    m_timer.setInterval(0);
    connect(&m_timer, &QTimer::timeout, this, &ProblemBadge::refresh);
    connect(&m_context, &forms::FormContext::problemsChanged, &m_timer, qOverload<>(&QTimer::start));
    hide();
    refresh();
}

void ProblemBadge::refresh() {
    m_timer.stop();
    const std::vector<core::Problem> mine = problemsOn(m_context.allProblems(), toStdString(m_page));
    ProblemCount tally;
    for (const core::Problem& p : mine) {
        (p.level == "warning" ? tally.warnings : tally.errors) += 1;
    }
    const int n = tally.total();
    setText(n < 100 ? QString::number(n) : QStringLiteral("99+"));
    setStyleSheet(QStringLiteral("QLabel { background-color: %1; color: #ffffff; border-radius: 8px; "
                                 "padding: 0px 4px; min-width: 8px; }")
                      .arg(QLatin1String(tally.errors ? kErrorBadge : kWarningBadge)));
    setToolTip(summaryText(mine));
    setAccessibleName(n ? QStringLiteral("%1: %2").arg(toQString(core::problemPageTitle(toStdString(m_page))),
                                                       countText(tally))
                        : QString());
    setProperty("level", n ? tally.level() : QString());
    setVisible(n > 0);
    const bool changed = tally.errors != m_count.errors || tally.warnings != m_count.warnings;
    m_count = tally;
    if (changed) {
        emit countChanged(tally.errors, tally.warnings);
    }
}

// ---- ProblemsPanel ---------------------------------------------------------------------------

ProblemsPanel::ProblemsPanel(forms::FormContext& context, QWidget* parent) : QWidget(parent), m_context(context) {
    setObjectName(QStringLiteral("problems"));
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_summary = new QLabel(this);
    m_summary->setObjectName(QStringLiteral("problems:summary"));
    layout->addWidget(m_summary);
    m_list = new QListWidget(this);
    m_list->setObjectName(QStringLiteral("problems:list"));
    m_list->setAccessibleName(QStringLiteral("Problems of the settings"));
    m_list->setToolTip(QStringLiteral("Problems of the settings: click one to go to its field"));
    m_list->setWordWrap(true);
    m_list->setTextElideMode(Qt::ElideNone);
    m_list->setSelectionMode(QAbstractItemView::SingleSelection);
    m_list->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
    m_list->setUniformItemSizes(false);
    layout->addWidget(m_list, 1);
    connect(m_list, &QListWidget::itemClicked, this, &ProblemsPanel::activate);
    connect(m_list, &QListWidget::itemActivated, this, &ProblemsPanel::activate);  // Enter, double-click
    m_timer.setSingleShot(true);
    m_timer.setInterval(0);
    connect(&m_timer, &QTimer::timeout, this, &ProblemsPanel::refresh);
    connect(&m_context, &forms::FormContext::problemsChanged, &m_timer, qOverload<>(&QTimer::start));
    refresh();
}

void ProblemsPanel::refresh() {
    m_timer.stop();
    std::vector<core::Problem> found = m_context.allProblems();
    const ProblemCount total = [&] {
        ProblemCount c;
        for (const core::Problem& p : found) {
            (p.level == "warning" ? c.warnings : c.errors) += 1;
        }
        return c;
    }();
    m_summary->setText(total.total() ? countText(total) : QStringLiteral("No problems"));
    if (found != m_problems || m_list->count() != static_cast<int>(found.size())) {
        m_problems = std::move(found);
        const int current = m_list->currentRow();
        m_list->clear();
        for (const core::Problem& p : m_problems) {
            const QString location = p.location.empty() ? QStringLiteral("Project") : toQString(p.location);
            const QString text = gui::displayPathsIn(toQString(p.text()));  // paths shown natively (WN-12)
            auto* item = new QListWidgetItem(levelIcon(this, p.level),
                                             location + QLatin1Char('\n') + gui::displayPathsIn(toQString(p.message)));
            item->setData(FieldRole, toQString(p.field));
            item->setData(PageRole, toQString(p.page()));
            item->setData(LevelRole, toQString(p.level));
            item->setToolTip(text);
            item->setData(Qt::AccessibleTextRole, QStringLiteral("%1: %2").arg(p.level == "warning" ? QStringLiteral("Warning")
                                                                                                    : QStringLiteral("Error"),
                                                                               text));
            m_list->addItem(item);
        }
        if (current >= 0 && current < m_list->count()) {
            m_list->setCurrentRow(current, QItemSelectionModel::NoUpdate);
        }
    }
    emit problemsShown(total.errors, total.warnings);
}

ProblemCount ProblemsPanel::count() const {
    ProblemCount c;
    for (const core::Problem& p : m_problems) {
        (p.level == "warning" ? c.warnings : c.errors) += 1;
    }
    return c;
}

ProblemCount ProblemsPanel::pageCount(const QString& page) const { return countOn(m_problems, toStdString(page)); }

QString ProblemsPanel::pageSummary(const QString& page) const {
    return summaryText(problemsOn(m_problems, toStdString(page)));
}

const core::Problem* ProblemsPanel::nextOnPage(const QString& page) {
    std::vector<const core::Problem*> mine;
    for (const core::Problem& p : m_problems) {
        if (toQString(p.page()) == page) {
            mine.push_back(&p);
        }
    }
    if (mine.empty()) {
        return nullptr;
    }
    int& k = m_next[page];
    const core::Problem* out = mine[static_cast<std::size_t>(k) % mine.size()];
    k = static_cast<int>((static_cast<std::size_t>(k) % mine.size()) + 1);
    return out;
}

bool ProblemsPanel::revealProblem(const QString& field) { return m_reveal && !field.isEmpty() && m_reveal(field); }

void ProblemsPanel::activate(QListWidgetItem* item) {
    if (!item) {
        return;
    }
    const QString field = item->data(FieldRole).toString();
    emit problemActivated(field, item->data(PageRole).toString());
    revealProblem(field);
}

}  // namespace mdbridge::app::problems
