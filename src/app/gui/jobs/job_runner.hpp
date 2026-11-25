#pragma once
// One background job at a time with a progress text, Cancel and the question when the window
// closes during a job (MDATK's gui/worker.py and 0.1.1's gui/jobs.py JobControl; AP-14).
//
// A job is `R job(const core::Reporter&)`; it runs on a worker thread of its own. Its log lines go
// to the runner's log function (called on the worker thread: it must be thread-safe, e.g.
// LogView::sink()), its progress reaches the interface thread throttled as MDATK's (a new fraction
// at least 0.5 % above the last one shown, 1.0, or a step back) and coalesced (at most one pending
// update; the newest wins). The result (`onSuccess(R)`) and failures are delivered on the thread
// the runner lives in. The texts are MDATK's: `<title> …`, `NN%  message`, `Cancelling …`,
// `Cancelled`, `Failed`, `<title>: done`; the log gets `--- <title> ---`, `Cancelled.` and
// `ERROR: <title> failed: <message>`; errorReported() is the window's error dialog (the line is
// already logged).

#include <chrono>
#include <exception>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include <QList>
#include <QObject>
#include <QPointer>
#include <QString>

#include "mdbridge/core/progress.hpp"

class QWidget;

namespace mdbridge::app::jobs {

enum class JobOutcome { Succeeded, Failed, Cancelled };

// A job that ended with an exception other than core::Cancelled.
struct JobFailure {
    QString title;              // the job's title
    QString message;            // what() of the exception ("unknown error" for a non-std one)
    std::exception_ptr error;   // the exception, for handlers that tell kinds apart
};

// The questions of close-while-running (AP-14) and their answers.
enum class CloseQuestion { JobRunning, JobNotStopped };
enum class CloseAnswer { Cancel, Close, KeepWaiting };

class JobRunner : public QObject {
    Q_OBJECT

public:
    using LogFunction = std::function<void(const QString& line)>;
    using FailureHandler = std::function<void(const JobFailure& failure)>;
    // Asks `question` about the job `title`; the default shows a QMessageBox over `parent`.
    using Asker = std::function<CloseAnswer(QWidget* parent, CloseQuestion question, const QString& title)>;

    // `log` receives every log line (the job's from its thread, the runner's own from this one).
    explicit JobRunner(LogFunction log, QObject* parent = nullptr);
    // A job still running is cancelled and waited for (unless the user chose "Close anyway").
    ~JobRunner() override;
    JobRunner(const JobRunner&) = delete;
    JobRunner& operator=(const JobRunner&) = delete;

    // Starts `job` unless one is running (then `Another job is still running.` is reported and
    // false returned). `onSuccess(result)` (or `onSuccess()` for a void job) and `onFailure` run on
    // this object's thread; onFailure runs before errorReported. The callables must be copyable.
    template <class Job, class Success>
    bool start(const QString& title, Job job, Success onSuccess, FailureHandler onFailure = {}) {
        using Result = std::invoke_result_t<Job&, const core::Reporter&>;
        Work work = [job = std::move(job), onSuccess = std::move(onSuccess)](const core::Reporter& reporter) mutable
            -> Delivery {
            if constexpr (std::is_void_v<Result>) {
                job(reporter);
                return [onSuccess]() mutable { onSuccess(); };
            } else {
                auto result = std::make_shared<Result>(job(reporter));
                return [onSuccess, result]() mutable { onSuccess(std::move(*result)); };
            }
        };
        return startWork(title, std::move(work), std::move(onFailure));
    }

    bool running() const { return m_run != nullptr; }
    QString title() const;   // of the running job ("" when none)
    double fraction() const { return m_fraction; }
    QString statusText() const { return m_text; }   // the progress bar's text ("Idle" at first)

    // Asks the running job to stop (`Cancelling …`); false when no job runs.
    bool cancel();
    // Waits up to `wait` for the running job to end, processing events (log lines, repaints; no
    // user input) meanwhile; true when no job runs any more.
    bool waitForFinished(std::chrono::milliseconds wait);
    // cancel() and waitForFinished(wait).
    bool stop(std::chrono::milliseconds wait);

    // Close-while-running (AP-14): true when the window may close. With a job running, asks first;
    // on "Close", cancels the job and waits for it (its own clean-up included) up to stopWait(),
    // then asks again whether to keep waiting or close anyway (the job is then abandoned).
    bool confirmClose(QWidget* parent);
    void setAsker(Asker asker) { m_asker = std::move(asker); }
    void setStopWait(std::chrono::milliseconds wait) { m_stopWait = wait; }
    std::chrono::milliseconds stopWait() const { return m_stopWait; }

    // Widgets or actions disabled while a job runs (Render, …) and enabled only while one runs
    // (Cancel, M7's Stop scripts), as 0.1.1's busy_widgets / stop_widgets.
    void addBusy(QObject* widgetOrAction);
    void addStoppable(QObject* widgetOrAction);

signals:
    void runningChanged(bool running);
    void statusChanged(double fraction, const QString& text);
    // The window's error dialog: `text` with `detail` (the line was logged already).
    void errorReported(const QString& text, const QString& detail);
    void finished(mdbridge::app::jobs::JobOutcome outcome);

private:
    using Delivery = std::function<void()>;
    using Work = std::function<Delivery(const core::Reporter&)>;
    struct Run;
    struct Channel;

    bool startWork(const QString& title, Work work, FailureHandler onFailure);
    void deliverProgress(Run* run);
    void finish(const std::shared_ptr<Run>& run);
    void setStatus(double fraction, const QString& text);
    void updateWidgets(bool busy);
    CloseAnswer ask(QWidget* parent, CloseQuestion question) const;

    LogFunction m_log;
    std::shared_ptr<Channel> m_channel;
    std::shared_ptr<Run> m_run;
    double m_fraction = 0.0;
    QString m_text;
    QList<QPointer<QObject>> m_busy;
    QList<QPointer<QObject>> m_stoppable;
    Asker m_asker;
    std::chrono::milliseconds m_stopWait{15000};   // 0.1.1 RNAdo: jobs.stop(wait=15.0)
    bool m_closing = false;
};

}  // namespace mdbridge::app::jobs
