#pragma once
// The log (MDATK's "Log (pipeline and cpptraj output)", app.py:369-382, with 0.1.1's gui/logpanel.py
// behaviour): lines arrive from any thread and are inserted in batches every 40 ms; a batch holds
// at most 2 000 lines, the newest, and a note counts the lines skipped when they come faster than
// that; the view keeps 5 000 to 6 000 lines (the oldest go first). Lines starting with ERROR are
// red and bold, WARNING / WAITING orange and bold, program stderr red, status lines bold, info
// lines blue (0.1.1's tags; lighter variants on a dark palette). Header: fold toggle (foldAction(),
// Ctrl+L, which the window adds to itself), title, line count, follow, Copy, Save, Clear.

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include <QString>
#include <QTextCharFormat>
#include <QTimer>
#include <QWidget>

class QAction;
class QLabel;
class QPlainTextEdit;

namespace mdbridge::app::jobs {

class LogView : public QWidget {
    Q_OBJECT

public:
    enum class Kind { Plain, Info, Stderr, Status };

    static constexpr int MaxLines = 6000;   // beyond this the oldest lines go ...
    static constexpr int KeepLines = 5000;  // ... down to this many
    static constexpr int FlushMs = 40;      // lines arriving within this time are inserted together
    static constexpr int MaxBatch = 2000;   // lines one flush inserts at most (the newest)

    explicit LogView(QWidget* parent = nullptr);
    ~LogView() override;
    LogView(const LogView&) = delete;
    LogView& operator=(const LogView&) = delete;

    // A message of one or more lines (a final newline ends the last line); safe from any thread.
    void feed(const QString& text, Kind kind = Kind::Plain);
    // feed() as a function that may outlive the view (later calls do nothing); for
    // JobRunner's log and core::Reporter callbacks.
    std::function<void(const QString&)> sink(Kind kind = Kind::Plain) const;

    void flush();                 // insert the queued lines now
    int lineCount() const { return m_lines; }
    QString text();               // the lines shown, each ended by a newline (flushes first)

    void clear();
    void copy();
    // Writes text() to `path` (asks for a name when empty); the path written, or "" on cancel /
    // failure (the log then says why).
    QString saveAs(const QString& path = QString());

    bool expanded() const;
    void setExpanded(bool expanded);
    bool following() const;
    void setFollowing(bool follow);

    QAction* foldAction() const { return m_fold; }   // checkable (checked = shown), Ctrl+L
    QAction* clearAction() const { return m_clear; }
    QAction* copyAction() const { return m_copy; }
    QAction* saveAction() const { return m_save; }
    QPlainTextEdit* textView() const { return m_view; }
    // Applied to each line as it is shown, on the interface thread (the window shows the paths in
    // lines natively with it, WN-12); none by default.
    void setLineTransform(std::function<QString(const QString&)> transform) { m_transform = std::move(transform); }

signals:
    void expandedChanged(bool expanded);

protected:
    void changeEvent(QEvent* event) override;

private:
    struct Queue;
    struct Line {
        QString text;
        Kind kind;
    };

    void showBatch(std::vector<Line> lines, long long skipped);
    const QTextCharFormat& formatFor(const Line& line) const;
    void updateFormats();
    void updateCount();
    void scrollToEnd();

    std::shared_ptr<Queue> m_queue;
    QTimer m_timer;
    QPlainTextEdit* m_view = nullptr;
    QLabel* m_count = nullptr;
    QAction* m_fold = nullptr;
    QAction* m_follow = nullptr;
    QAction* m_clear = nullptr;
    QAction* m_copy = nullptr;
    QAction* m_save = nullptr;
    std::function<QString(const QString&)> m_transform;
    int m_lines = 0;
    QTextCharFormat m_plain, m_info, m_stderr, m_status, m_error, m_waiting;

    friend struct Queue;
};

}  // namespace mdbridge::app::jobs
