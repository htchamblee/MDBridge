#include "gui/jobs/job_runner.hpp"

#include <algorithm>
#include <mutex>
#include <string>
#include <thread>

#include <QAbstractButton>
#include <QAction>
#include <QEventLoop>
#include <QMessageBox>
#include <QMetaObject>
#include <QPushButton>
#include <QTimer>
#include <QWidget>

#include "mdbridge/core/errors.hpp"

namespace mdbridge::app::jobs {

// The runner as seen from worker threads: cleared when the runner goes, so a job that outlives it
// (abandoned at close) posts into the void.
struct JobRunner::Channel {
    std::mutex mutex;
    JobRunner* receiver = nullptr;

    template <class Fn>
    void post(Fn fn) {
        const std::lock_guard<std::mutex> lock(mutex);
        if (receiver) {
            // the receiver's posted events are discarded if it is deleted before they run
            QMetaObject::invokeMethod(receiver, [r = receiver, fn = std::move(fn)]() mutable { fn(r); },
                                      Qt::QueuedConnection);
        }
    }
};

struct JobRunner::Run {
    QString title;
    FailureHandler onFailure;
    core::Reporter reporter;
    std::thread thread;
    bool abandoned = false;  // "Close anyway": detached, never joined

    // progress throttling (worker threads) and the pending update (both threads)
    std::mutex mutex;
    double last = -1.0;
    bool progressPending = false;
    double pendingFraction = 0.0;
    std::string pendingMessage;

    // the outcome, written by the worker before it posts finish()
    JobOutcome outcome = JobOutcome::Failed;
    JobFailure failure;
    Delivery delivery;
};

JobRunner::JobRunner(LogFunction log, QObject* parent)
    : QObject(parent), m_log(std::move(log)), m_channel(std::make_shared<Channel>()), m_text(tr("Idle")) {
    m_channel->receiver = this;
}

JobRunner::~JobRunner() {
    {
        const std::lock_guard<std::mutex> lock(m_channel->mutex);
        m_channel->receiver = nullptr;
    }
    if (m_run) {
        m_run->reporter.cancel();
        if (m_run->abandoned) {
            m_run->thread.detach();
        } else if (m_run->thread.joinable()) {
            m_run->thread.join();  // AP-14: never end the process under a job's feet
        }
    }
}

QString JobRunner::title() const { return m_run ? m_run->title : QString(); }

bool JobRunner::startWork(const QString& title, Work work, FailureHandler onFailure) {
    if (m_run) {
        const QString text = tr("Another job is still running.");
        if (m_log) m_log(QStringLiteral("ERROR: ") + text);
        emit errorReported(text, QString());
        return false;
    }
    auto run = std::make_shared<Run>();
    run->title = title;
    run->onFailure = std::move(onFailure);
    Run* raw = run.get();
    const std::shared_ptr<Channel> channel = m_channel;
    core::Reporter::LogCallback logLine;
    if (m_log) {
        logLine = [log = m_log](const std::string& message) { log(QString::fromStdString(message)); };
    }
    run->reporter = core::Reporter(std::move(logLine), [raw, channel](double fraction, const std::string& message) {
        const std::lock_guard<std::mutex> lock(raw->mutex);
        // MDATK's worker.py: a new value at least 0.5 % on, the end, or a step back
        if (!(fraction >= 1.0 || fraction - raw->last >= 0.005 || fraction < raw->last)) return;
        raw->last = fraction;
        raw->pendingFraction = fraction;
        raw->pendingMessage = message;
        if (raw->progressPending) return;  // one update in flight at a time: it takes the newest
        raw->progressPending = true;
        channel->post([raw](JobRunner* receiver) { receiver->deliverProgress(raw); });
    });

    m_run = run;
    setStatus(0.0, tr("%1 …").arg(title));
    updateWidgets(true);
    if (m_log) m_log(QStringLiteral("--- %1 ---").arg(title));
    emit runningChanged(true);

    run->thread = std::thread([run, channel, work = std::move(work)]() mutable {
        try {
            run->reporter.checkCancel();  // cancelled before it began
            run->delivery = work(run->reporter);
            run->outcome = JobOutcome::Succeeded;
        } catch (const core::Cancelled&) {
            run->outcome = JobOutcome::Cancelled;
        } catch (const std::exception& error) {
            run->outcome = JobOutcome::Failed;
            run->failure = {run->title, QString::fromStdString(error.what()), std::current_exception()};
        } catch (...) {
            run->outcome = JobOutcome::Failed;
            run->failure = {run->title, QStringLiteral("unknown error"), std::current_exception()};
        }
        work = {};  // the job's captures go here, on its own thread
        channel->post([run](JobRunner* receiver) { receiver->finish(run); });
    });
    return true;
}

void JobRunner::deliverProgress(Run* run) {
    if (!m_run || m_run.get() != run) return;  // a late update of a finished job
    double fraction = 0.0;
    std::string message;
    {
        const std::lock_guard<std::mutex> lock(run->mutex);
        run->progressPending = false;
        fraction = run->pendingFraction;
        message = run->pendingMessage;
    }
    fraction = std::clamp(fraction, 0.0, 1.0);
    QString text = m_text;
    if (!message.empty()) {
        text = QString::asprintf("%4.0f%%  ", fraction * 100.0) + QString::fromStdString(message);
    }
    setStatus(fraction, text);
}

void JobRunner::finish(const std::shared_ptr<Run>& run) {
    if (run->thread.joinable()) run->thread.join();  // it has posted this as its last act
    if (run != m_run) return;
    m_run.reset();
    updateWidgets(false);
    emit runningChanged(false);

    JobOutcome outcome = run->outcome;
    JobFailure failure = run->failure;
    if (outcome == JobOutcome::Succeeded) {
        setStatus(1.0, tr("%1: done").arg(run->title));
        try {
            if (run->delivery) run->delivery();
        } catch (const std::exception& error) {  // the result handler failed: the job did
            outcome = JobOutcome::Failed;
            failure = {run->title, QString::fromStdString(error.what()), std::current_exception()};
        } catch (...) {
            outcome = JobOutcome::Failed;
            failure = {run->title, QStringLiteral("unknown error"), std::current_exception()};
        }
    }
    if (outcome == JobOutcome::Cancelled) {
        setStatus(m_fraction, tr("Cancelled"));
        if (m_log) m_log(tr("Cancelled."));
    } else if (outcome == JobOutcome::Failed) {
        setStatus(m_fraction, tr("Failed"));
        if (run->onFailure) run->onFailure(failure);
        const QString text = tr("%1 failed").arg(run->title);
        if (m_log) m_log(QStringLiteral("ERROR: %1: %2").arg(text, failure.message));
        emit errorReported(text, failure.message);
    }
    emit finished(outcome);
}

bool JobRunner::cancel() {
    if (!m_run) return false;
    m_run->reporter.cancel();
    setStatus(m_fraction, tr("Cancelling …"));
    return true;
}

bool JobRunner::waitForFinished(std::chrono::milliseconds wait) {
    if (!m_run) return true;
    QEventLoop loop;
    const auto connection = connect(this, &JobRunner::finished, &loop, &QEventLoop::quit);
    QTimer::singleShot(wait, &loop, &QEventLoop::quit);
    loop.exec(QEventLoop::ExcludeUserInputEvents);
    disconnect(connection);
    return !m_run;
}

bool JobRunner::stop(std::chrono::milliseconds wait) {
    cancel();
    return waitForFinished(wait);
}

bool JobRunner::confirmClose(QWidget* parent) {
    if (!m_run) return true;
    if (m_closing) return false;  // a second close request while the first one waits
    m_closing = true;
    struct Reset {
        bool& flag;
        ~Reset() { flag = false; }
    } reset{m_closing};

    if (ask(parent, CloseQuestion::JobRunning) != CloseAnswer::Close) return false;
    if (!m_run) return true;
    if (stop(m_stopWait)) return true;
    for (;;) {
        const CloseAnswer answer = ask(parent, CloseQuestion::JobNotStopped);
        if (!m_run) return true;  // it stopped while the question was shown
        if (answer == CloseAnswer::KeepWaiting) {
            if (waitForFinished(m_stopWait)) return true;
        } else if (answer == CloseAnswer::Close) {
            m_run->abandoned = true;
            if (m_log) m_log(tr("WARNING: closed while %1 had not stopped").arg(m_run->title));
            return true;
        } else {
            return false;  // the window stays; the job goes on stopping
        }
    }
}

CloseAnswer JobRunner::ask(QWidget* parent, CloseQuestion question) const {
    const QString jobTitle = title();
    if (m_asker) return m_asker(parent, question, jobTitle);
    QMessageBox box(parent);
    box.setIcon(QMessageBox::Warning);
    QAbstractButton* close = nullptr;
    QAbstractButton* keep = nullptr;
    if (question == CloseQuestion::JobRunning) {
        box.setWindowTitle(tr("A job is still running"));
        box.setText(tr("%1 is still running.").arg(jobTitle));
        box.setInformativeText(tr("Closing cancels it and waits until it has stopped."));
        box.addButton(QMessageBox::Cancel);
        close = box.addButton(tr("Cancel it and close"), QMessageBox::AcceptRole);
        box.setDefaultButton(QMessageBox::Cancel);
    } else {
        box.setWindowTitle(tr("The job has not stopped yet"));
        box.setText(tr("%1 has not stopped yet.").arg(jobTitle));
        box.setInformativeText(tr("It may still be finishing or removing the files it writes. Closing now "
                                  "can leave incomplete files behind."));
        keep = box.addButton(tr("Keep waiting"), QMessageBox::AcceptRole);
        close = box.addButton(tr("Close anyway"), QMessageBox::DestructiveRole);
        box.addButton(QMessageBox::Cancel);
        box.setDefaultButton(qobject_cast<QPushButton*>(keep));
    }
    box.exec();
    if (box.clickedButton() == close) return CloseAnswer::Close;
    if (keep != nullptr && box.clickedButton() == keep) return CloseAnswer::KeepWaiting;
    return CloseAnswer::Cancel;
}

void JobRunner::setStatus(double fraction, const QString& text) {
    m_fraction = fraction;
    m_text = text;
    emit statusChanged(fraction, text);
}

namespace {

void setEnabled(QObject* object, bool enabled) {
    if (auto* widget = qobject_cast<QWidget*>(object)) {
        widget->setEnabled(enabled);
    } else if (auto* action = qobject_cast<QAction*>(object)) {
        action->setEnabled(enabled);
    }
}

}  // namespace

void JobRunner::updateWidgets(bool busy) {
    for (const QPointer<QObject>& object : std::as_const(m_busy)) {
        if (object) setEnabled(object, !busy);
    }
    for (const QPointer<QObject>& object : std::as_const(m_stoppable)) {
        if (object) setEnabled(object, busy);
    }
}

void JobRunner::addBusy(QObject* widgetOrAction) {
    m_busy.append(widgetOrAction);
    setEnabled(widgetOrAction, !running());
}

void JobRunner::addStoppable(QObject* widgetOrAction) {
    m_stoppable.append(widgetOrAction);
    setEnabled(widgetOrAction, running());
}

}  // namespace mdbridge::app::jobs
