#include "gui/tools/tool_window.hpp"

#include <algorithm>
#include <utility>

#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QStyle>
#include <QTextOption>
#include <QToolButton>
#include <QVBoxLayout>

#include "gui/forms/field_widgets.hpp"
#include "gui/tools/host.hpp"
#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::app::tools {

namespace {

constexpr int kFormWidth = 500;      // MDATK: set_size_request(500, -1)
constexpr int kSplitterAt = 520;     // MDATK: paned.set_position(520)
constexpr int kComboChars = 16;      // long labels are elided, not widening the form

}  // namespace

// ---- PathEntry -----------------------------------------------------------------------------

PathEntry::PathEntry(const QString& name, const QString& dialogTitle, forms::PathMode mode, QStringList filters,
                     QWidget* parent)
    : QWidget(parent), m_title(dialogTitle), m_mode(mode), m_filters(std::move(filters)) {
    setObjectName(name + QStringLiteral(":field"));
    auto* row = new QHBoxLayout(this);
    row->setContentsMargins(0, 0, 0, 0);
    row->setSpacing(4);
    m_edit = new QLineEdit(this);
    m_edit->setObjectName(name);
    m_edit->setMinimumWidth(120);
    row->addWidget(m_edit, 1);
    m_browse = new QToolButton(this);
    m_browse->setObjectName(name + QStringLiteral(":browse"));
    m_browse->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    m_browse->setText(QStringLiteral("Browse…"));
    m_browse->setToolTip(QStringLiteral("Browse: ") + dialogTitle);
    m_browse->setAccessibleName(QStringLiteral("Browse: ") + dialogTitle);
    m_browse->setFocusPolicy(Qt::StrongFocus);
    row->addWidget(m_browse);
    connect(m_browse, &QToolButton::clicked, this, &PathEntry::browse);
    connect(m_edit, &QLineEdit::textChanged, this, &PathEntry::pathChanged);
}

QString PathEntry::path() const { return m_edit->text().trimmed(); }

void PathEntry::setPath(const QString& text) { m_edit->setText(text); }

void PathEntry::browse() {
    forms::FileDialogRequest request;
    request.mode = m_mode;
    request.title = m_title;
    const QString current = path();
    request.start = current.isEmpty() ? QDir::currentPath() : QFileInfo(current).absoluteFilePath();
    request.filters = m_filters;
    request.confirmOverwrite = true;
    const QStringList picked = forms::runFileDialog(window(), request);
    if (!picked.isEmpty()) {
        m_edit->setText(QDir::toNativeSeparators(picked.front()));  // as MDATK: the file chosen, absolute
    }
}

// ---- ToolWindow ----------------------------------------------------------------------------

ToolWindow::ToolWindow(Host& host, QString id, const QString& title, QWidget* parent)
    : QWidget(parent, Qt::Window), m_host(host), m_id(std::move(id)) {
    setObjectName(m_id);
    setWindowTitle(title);
    setAttribute(Qt::WA_DeleteOnClose);
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(0, 0, 0, 0);
    outer->setSpacing(0);
    m_splitter = new QSplitter(Qt::Horizontal, this);
    m_splitter->setObjectName(nameOf(QStringLiteral("splitter")));
    m_splitter->setChildrenCollapsible(false);
    outer->addWidget(m_splitter, 1);
    auto* bar = new QWidget(this);
    bar->setObjectName(nameOf(QStringLiteral("bar")));
    m_bar = new QHBoxLayout(bar);
    m_bar->setContentsMargins(8, 8, 8, 8);
    m_bar->setSpacing(6);
    m_bar->addStretch(1);
    outer->addWidget(bar);
}

ToolWindow::~ToolWindow() = default;

QWidget* ToolWindow::control(const QString& key) const {
    const QList<QWidget*> found = findChildren<QWidget*>(nameOf(key));
    return found.size() == 1 ? found.front() : nullptr;
}

void ToolWindow::showReport(const QString& text) { m_report->setPlainText(text); }

QScrollArea* ToolWindow::makeFormScroll(const QString& key, const QString& accessibleName, QVBoxLayout** column,
                                        QWidget* parent) {
    auto* scroll = new QScrollArea(parent);
    scroll->setObjectName(nameOf(key) + QStringLiteral(":scroll"));
    scroll->setWidgetResizable(true);
    scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setAccessibleName(accessibleName);  // it takes the focus to scroll the form with the keys
    scroll->setMinimumWidth(kFormWidth);
    auto* body = new QWidget;
    body->setObjectName(nameOf(key) + QStringLiteral(":body"));
    auto* layout = new QVBoxLayout(body);
    layout->setContentsMargins(10, 10, 10, 10);
    layout->setSpacing(8);
    scroll->setWidget(body);
    *column = layout;
    return scroll;
}

QVBoxLayout* ToolWindow::addReport() {
    auto* right = new QWidget(m_splitter);
    right->setObjectName(nameOf(QStringLiteral("report_pane")));
    auto* layout = new QVBoxLayout(right);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(4);
    auto* heading = new QLabel(QStringLiteral("<b>Report</b>"), right);
    heading->setObjectName(nameOf(QStringLiteral("report_heading")));
    layout->addWidget(heading);
    m_report = new QPlainTextEdit(right);
    m_report->setObjectName(nameOf(QStringLiteral("report")));
    m_report->setReadOnly(true);
    m_report->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    m_report->setFont(render::interfaceFixedFont(QApplication::font().pointSizeF()));
    m_report->setWordWrapMode(QTextOption::WrapAtWordBoundaryOrAnywhere);  // GTK's WORD_CHAR
    m_report->setAccessibleName(QStringLiteral("Report"));
    m_report->setTabChangesFocus(true);
    heading->setBuddy(m_report);
    layout->addWidget(m_report, 1);
    m_splitter->addWidget(right);
    return layout;
}

QPushButton* ToolWindow::addButton(const QString& key, const QString& text, const QString& tooltip) {
    auto* button = new QPushButton(text, m_bar->parentWidget());
    button->setObjectName(nameOf(key));
    if (!tooltip.isEmpty()) button->setToolTip(tooltip);
    button->setAutoDefault(false);
    m_bar->addWidget(button);
    return button;
}

void ToolWindow::finishLayout(int width, int height) {
    m_splitter->setStretchFactor(0, 0);
    m_splitter->setStretchFactor(1, 1);
    resize(width, height);
    m_splitter->setSizes({kSplitterAt, width - kSplitterAt});
}

QLineEdit* ToolWindow::makeLine(const QString& key, const QString& accessibleName, const QString& placeholder,
                                const QString& tooltip, int widthChars) {
    auto* edit = new QLineEdit(this);
    edit->setObjectName(nameOf(key));
    edit->setAccessibleName(accessibleName);
    if (!placeholder.isEmpty()) edit->setPlaceholderText(placeholder);
    if (!tooltip.isEmpty()) edit->setToolTip(tooltip);
    if (widthChars > 0) {  // MDATK's W.entry(width=...): a short field (chain ids)
        edit->setFixedWidth(edit->fontMetrics().horizontalAdvance(QString(widthChars + 2, QLatin1Char('M'))));
    }
    return edit;
}

QComboBox* ToolWindow::makeCombo(const QString& key, const QString& accessibleName,
                                 std::initializer_list<std::pair<const char*, const char*>> items, const char* active) {
    auto* combo = new QComboBox(this);
    combo->setObjectName(nameOf(key));
    combo->setAccessibleName(accessibleName);
    combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    combo->setMinimumContentsLength(kComboChars);
    for (const auto& [value, label] : items) {
        const QString text = QString::fromUtf8(label);
        combo->addItem(text, QString::fromUtf8(value));
        combo->setItemData(combo->count() - 1, text, Qt::ToolTipRole);  // the whole label of an elided item
    }
    combo->setCurrentIndex(std::max(0, combo->findData(QString::fromUtf8(active))));
    forms::guardWheel(combo);
    return combo;
}

QCheckBox* ToolWindow::makeCheck(const QString& key, const QString& text, bool checked, const QString& tooltip) {
    auto* check = new QCheckBox(text, this);
    check->setObjectName(nameOf(key));
    check->setChecked(checked);
    if (!tooltip.isEmpty()) check->setToolTip(tooltip);
    return check;
}

QSpinBox* ToolWindow::makeSpin(const QString& key, const QString& accessibleName, int lo, int hi, int value) {
    auto* spin = new QSpinBox(this);
    spin->setObjectName(nameOf(key));
    spin->setAccessibleName(accessibleName);
    spin->setRange(lo, hi);
    spin->setValue(value);
    spin->setMinimumWidth(84);
    forms::guardWheel(spin);
    return spin;
}

QDoubleSpinBox* ToolWindow::makeDoubleSpin(const QString& key, const QString& accessibleName, double lo, double hi,
                                           double step, int decimals, double value) {
    auto* spin = new QDoubleSpinBox(this);
    spin->setObjectName(nameOf(key));
    spin->setAccessibleName(accessibleName);
    spin->setDecimals(decimals);
    spin->setRange(lo, hi);
    spin->setSingleStep(step);
    spin->setValue(value);
    spin->setMinimumWidth(84);
    forms::guardWheel(spin);
    return spin;
}

PathEntry* ToolWindow::makePath(const QString& key, const QString& accessibleName, const QString& dialogTitle,
                                forms::PathMode mode, const std::vector<forms::FileFilter>& filters) {
    QStringList names;
    for (const forms::FileFilter& f : filters) names << forms::nameFilter(f);
    if (!names.isEmpty()) names << QStringLiteral("All files (*)");
    auto* entry = new PathEntry(nameOf(key), dialogTitle, mode, names, this);
    entry->edit()->setAccessibleName(accessibleName);
    return entry;
}

QString ToolWindow::currentId(const QComboBox* combo) { return combo->currentData().toString(); }

void ToolWindow::closeEvent(QCloseEvent* event) {
    if (!mayClose()) {
        event->ignore();
        return;
    }
    QWidget::closeEvent(event);
}

}  // namespace mdbridge::app::tools
