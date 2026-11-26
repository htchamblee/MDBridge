#include "gui/jobs/log_view.hpp"

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <utility>
#include <vector>

#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QEvent>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QKeySequence>
#include <QLabel>
#include <QMetaObject>
#include <QPalette>
#include <QPlainTextEdit>
#include <QSaveFile>
#include <QScrollBar>
#include <QStringList>
#include <QTextCursor>
#include <QTextDocument>
#include <QThread>
#include <QToolButton>
#include <QVBoxLayout>

#include "mdbridge/render/text_qt.hpp"

namespace mdbridge::app::jobs {

// Lines waiting for the next flush, shared with the threads that feed them.
struct LogView::Queue {
    std::mutex mutex;
    LogView* view = nullptr;  // cleared when the view goes
    std::vector<Line> pending;
    long long skipped = 0;  // queued lines dropped since the last flush
    bool scheduled = false;

    void push(const QString& text, Kind kind) {
        QStringList parts = text.split(QLatin1Char('\n'));
        if (parts.size() > 1 && parts.back().isEmpty()) parts.removeLast();  // "line\n" is one line
        const std::lock_guard<std::mutex> lock(mutex);
        if (!view) return;
        for (QString& part : parts) {
            if (part.endsWith(QLatin1Char('\r'))) part.chop(1);
            pending.push_back({std::move(part), kind});
        }
        if (pending.size() > static_cast<std::size_t>(2 * MaxBatch)) dropOldest(pending.size() - static_cast<std::size_t>(MaxBatch));
        if (!scheduled) {
            scheduled = true;
            LogView* target = view;
            if (QThread::currentThread() == target->thread()) {
                target->m_timer.start();
            } else {
                QMetaObject::invokeMethod(target, [target] { target->m_timer.start(); }, Qt::QueuedConnection);
            }
        }
    }

    // in chunks, so the work per line stays constant however fast lines come (the lock is held)
    void dropOldest(std::size_t n) {
        pending.erase(pending.begin(), pending.begin() + static_cast<std::ptrdiff_t>(n));
        skipped += static_cast<long long>(n);
    }
};

LogView::LogView(QWidget* parent) : QWidget(parent), m_queue(std::make_shared<Queue>()) {
    m_queue->view = this;
    m_timer.setSingleShot(true);
    m_timer.setInterval(FlushMs);
    connect(&m_timer, &QTimer::timeout, this, &LogView::flush);

    m_fold = new QAction(tr("Show the log"), this);
    m_fold->setObjectName(QStringLiteral("log.fold"));
    m_fold->setCheckable(true);
    m_fold->setChecked(true);
    m_fold->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_L));
    m_fold->setShortcutContext(Qt::WindowShortcut);
    m_fold->setToolTip(tr("Show or hide the log (Ctrl+L)"));
    m_follow = new QAction(tr("Follow"), this);
    m_follow->setObjectName(QStringLiteral("log.follow"));
    m_follow->setCheckable(true);
    m_follow->setChecked(true);
    m_follow->setToolTip(tr("Follow new lines (auto-scroll)"));
    m_copy = new QAction(tr("Copy"), this);
    m_copy->setObjectName(QStringLiteral("log.copy"));
    m_copy->setToolTip(tr("Copy the lines shown"));
    m_save = new QAction(tr("Save…"), this);
    m_save->setObjectName(QStringLiteral("log.save"));
    m_save->setToolTip(tr("Save the lines shown as a text file"));
    m_clear = new QAction(tr("Clear"), this);
    m_clear->setObjectName(QStringLiteral("log.clear"));
    m_clear->setToolTip(tr("Clear the log"));

    auto* header = new QHBoxLayout;
    header->setContentsMargins(4, 4, 4, 4);
    header->setSpacing(6);
    auto* foldButton = new QToolButton(this);
    foldButton->setDefaultAction(m_fold);
    foldButton->setArrowType(Qt::DownArrow);
    foldButton->setAutoRaise(true);
    foldButton->setAccessibleName(tr("Show or hide the log"));
    header->addWidget(foldButton);
    auto* title = new QLabel(tr("<b>Log</b> (pipeline and cpptraj output)"), this);
    header->addWidget(title);
    m_count = new QLabel(this);
    m_count->setObjectName(QStringLiteral("log.count"));
    m_count->setEnabled(false);  // a hint
    header->addWidget(m_count, 1);
    for (QAction* action : {m_follow, m_copy, m_save, m_clear}) {
        auto* button = new QToolButton(this);
        button->setDefaultAction(action);
        button->setToolButtonStyle(Qt::ToolButtonTextOnly);
        header->addWidget(button);
    }

    m_view = new QPlainTextEdit(this);
    m_view->setObjectName(QStringLiteral("log.text"));
    m_view->setReadOnly(true);
    m_view->setUndoRedoEnabled(false);
    m_view->setLineWrapMode(QPlainTextEdit::WidgetWidth);
    m_view->setWordWrapMode(QTextOption::WrapAtWordBoundaryOrAnywhere);
    m_view->setFont(render::interfaceFixedFont(QApplication::font().pointSizeF()));
    m_view->setMinimumHeight(60);
    m_view->setAccessibleName(tr("Log"));

    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addLayout(header);
    layout->addWidget(m_view, 1);

    connect(m_fold, &QAction::toggled, this, [this, foldButton](bool on) {
        foldButton->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
        m_view->setVisible(on);
        emit expandedChanged(on);
    });
    connect(m_follow, &QAction::toggled, this, [this](bool on) {
        if (on) scrollToEnd();
    });
    connect(m_copy, &QAction::triggered, this, &LogView::copy);
    connect(m_save, &QAction::triggered, this, [this] { saveAs(); });
    connect(m_clear, &QAction::triggered, this, &LogView::clear);
    updateFormats();
    updateCount();
}

LogView::~LogView() {
    const std::lock_guard<std::mutex> lock(m_queue->mutex);
    m_queue->view = nullptr;
    m_queue->pending.clear();
}

void LogView::feed(const QString& text, Kind kind) { m_queue->push(text, kind); }

std::function<void(const QString&)> LogView::sink(Kind kind) const {
    return [queue = m_queue, kind](const QString& text) { queue->push(text, kind); };
}

void LogView::flush() {
    std::vector<Line> lines;
    long long skipped = 0;
    {
        const std::lock_guard<std::mutex> lock(m_queue->mutex);
        if (m_queue->pending.size() > static_cast<std::size_t>(MaxBatch)) {
            m_queue->dropOldest(m_queue->pending.size() - static_cast<std::size_t>(MaxBatch));
        }
        lines.swap(m_queue->pending);
        skipped = m_queue->skipped;
        m_queue->skipped = 0;
        m_queue->scheduled = false;
    }
    m_timer.stop();
    if (!lines.empty() || skipped > 0) showBatch(std::move(lines), skipped);
}

void LogView::showBatch(std::vector<Line> lines, long long skipped) {
    if (skipped > 0) {
        lines.insert(lines.begin(),
                     Line{tr("… %1 line(s) not shown: they came faster than the log can show them").arg(skipped),
                          Kind::Info});
    }
    QTextDocument* document = m_view->document();
    const int total = m_lines + static_cast<int>(lines.size());
    if (total > MaxLines) {  // the oldest lines go first: the view never holds more than MaxLines
        const int cut = total - KeepLines;
        const int old = std::min(cut, m_lines);
        if (old >= m_lines) {
            document->clear();
            m_lines = 0;
        } else if (old > 0) {
            QTextCursor cursor(document);
            cursor.movePosition(QTextCursor::Start);
            cursor.movePosition(QTextCursor::NextBlock, QTextCursor::KeepAnchor, old);
            cursor.removeSelectedText();
            m_lines -= old;
        }
        lines.erase(lines.begin(), lines.begin() + (cut - old));
    }
    if (m_transform) {
        for (Line& line : lines) line.text = m_transform(line.text);
    }
    if (!lines.empty()) {
        QTextCursor cursor(document);
        cursor.movePosition(QTextCursor::End);
        cursor.beginEditBlock();
        // runs of lines of one style go in with one call
        std::size_t i = 0;
        bool first = m_lines == 0;
        while (i < lines.size()) {
            const QTextCharFormat& format = formatFor(lines[i]);
            QString run;
            std::size_t j = i;
            for (; j < lines.size() && &formatFor(lines[j]) == &format; ++j) {
                if (!first) run += QLatin1Char('\n');
                first = false;
                run += lines[j].text;
            }
            cursor.insertText(run, format);
            i = j;
        }
        cursor.endEditBlock();
        m_lines += static_cast<int>(lines.size());
    }
    updateCount();
    if (m_follow->isChecked()) scrollToEnd();
}

const QTextCharFormat& LogView::formatFor(const Line& line) const {
    // 0.1.1 logpanel.py Line.tag()
    if (line.text.startsWith(QLatin1String("ERROR"))) return m_error;
    if (line.text.startsWith(QLatin1String("WAITING")) || line.text.startsWith(QLatin1String("WARNING"))) {
        return m_waiting;
    }
    switch (line.kind) {
    case Kind::Stderr: return m_stderr;
    case Kind::Status: return m_status;
    case Kind::Info: return m_info;
    case Kind::Plain: break;
    }
    return m_plain;
}

void LogView::updateFormats() {
    // 0.1.1's colours (GNOME palette) on a light base; lighter shades of the same hues on a dark one
    const bool dark = palette().color(QPalette::Base).lightness() < 128;
    m_plain = QTextCharFormat();
    m_info = QTextCharFormat();
    m_info.setForeground(QColor(dark ? QStringLiteral("#99c1f1") : QStringLiteral("#1a5fb4")));
    m_stderr = QTextCharFormat();
    m_stderr.setForeground(QColor(dark ? QStringLiteral("#f66151") : QStringLiteral("#a51d2d")));
    m_status = QTextCharFormat();
    m_status.setFontWeight(QFont::Bold);
    m_error = QTextCharFormat();
    m_error.setForeground(QColor(dark ? QStringLiteral("#ff7b63") : QStringLiteral("#c01c28")));
    m_error.setFontWeight(QFont::Bold);
    m_waiting = QTextCharFormat();
    m_waiting.setForeground(QColor(dark ? QStringLiteral("#ffa348") : QStringLiteral("#b56a00")));
    m_waiting.setFontWeight(QFont::Bold);
}

void LogView::changeEvent(QEvent* event) {
    if (event->type() == QEvent::PaletteChange) updateFormats();
    QWidget::changeEvent(event);
}

void LogView::updateCount() {
    m_count->setText(m_lines > 0 ? tr("%1 line(s)").arg(m_lines) : QString());
}

void LogView::scrollToEnd() {
    QScrollBar* bar = m_view->verticalScrollBar();
    bar->setValue(bar->maximum());
}

QString LogView::text() {
    flush();
    if (m_lines == 0) return QString();
    return m_view->document()->toPlainText() + QLatin1Char('\n');
}

void LogView::clear() {
    {
        const std::lock_guard<std::mutex> lock(m_queue->mutex);
        m_queue->pending.clear();
        m_queue->skipped = 0;
    }
    m_view->document()->clear();
    m_lines = 0;
    updateCount();
}

void LogView::copy() { QApplication::clipboard()->setText(text()); }

QString LogView::saveAs(const QString& path) {
    const QString content = text();
    QString target = path;
    if (target.isEmpty()) {
        target = QFileDialog::getSaveFileName(this, tr("Save the log"), QStringLiteral("mdbridge_log.txt"),
                                              tr("Text files (*.txt *.log)"));
        if (target.isEmpty()) return QString();
    }
    QSaveFile file(target);
    if (!file.open(QIODevice::WriteOnly) || file.write(content.toUtf8()) < 0 || !file.commit()) {
        feed(QStringLiteral("ERROR: ") + tr("could not write %1: %2").arg(target, file.errorString()));
        return QString();
    }
    feed(tr("Log saved to %1").arg(target));
    return target;
}

bool LogView::expanded() const { return m_fold->isChecked(); }
void LogView::setExpanded(bool expanded) { m_fold->setChecked(expanded); }
bool LogView::following() const { return m_follow->isChecked(); }
void LogView::setFollowing(bool follow) { m_follow->setChecked(follow); }

}  // namespace mdbridge::app::jobs
