// The job runner and the log view (src/app/gui/jobs): MDATK's job texts and one-job rule, the
// progress throttle, cancel before / during / after a job, failures delivered to the interface
// thread, close-while-running (AP-14), and a log that stays bounded and keeps the event loop
// responsive under 10^6 lines (offscreen platform).

#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QPalette>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollBar>
#include <QSignalSpy>
#include <QStringList>
#include <QTemporaryDir>
#include <QTest>
#include <QTextBlock>
#include <QTextCursor>
#include <QTextDocument>
#include <QTimer>

#include "gui/jobs/job_runner.hpp"
#include "gui/jobs/log_view.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/progress.hpp"

using mdbridge::app::jobs::CloseAnswer;
using mdbridge::app::jobs::CloseQuestion;
using mdbridge::app::jobs::JobFailure;
using mdbridge::app::jobs::JobOutcome;
using mdbridge::app::jobs::JobRunner;
using mdbridge::app::jobs::LogView;
using mdbridge::core::Reporter;
using namespace std::chrono_literals;

namespace {

// Log lines collected from any thread.
struct Lines {
    std::shared_ptr<std::mutex> mutex = std::make_shared<std::mutex>();
    std::shared_ptr<QStringList> lines = std::make_shared<QStringList>();
    JobRunner::LogFunction function() const {
        return [m = mutex, l = lines](const QString& line) {
            const std::lock_guard<std::mutex> lock(*m);
            l->append(line);
        };
    }
    QStringList get() const {
        const std::lock_guard<std::mutex> lock(*mutex);
        return *lines;
    }
};

// The longest time the event loop went without running a 5 ms timer while `busy()` holds.
template <class Busy>
double longestStallSeconds(Busy busy, std::chrono::seconds limit = 120s) {
    QElapsedTimer clock;
    clock.start();
    qint64 last = 0;
    qint64 longest = 0;
    QTimer beat;
    beat.setInterval(5);
    QObject::connect(&beat, &QTimer::timeout, [&] {
        const qint64 now = clock.nsecsElapsed();
        longest = std::max(longest, now - last);
        last = now;
    });
    beat.start();
    while (busy() && clock.elapsed() < std::chrono::milliseconds(limit).count()) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 5);
    }
    longest = std::max(longest, clock.nsecsElapsed() - last);
    return static_cast<double>(longest) / 1e9;
}

QStringList shownLines(LogView& log) {
    QString text = log.text();
    if (text.endsWith(QLatin1Char('\n'))) text.chop(1);
    return text.isEmpty() ? QStringList() : text.split(QLatin1Char('\n'));
}

QTextCharFormat formatOfLine(LogView& log, int line) {
    log.flush();
    const QTextBlock block = log.textView()->document()->findBlockByNumber(line);
    QTextCursor cursor(block);
    cursor.movePosition(QTextCursor::NextCharacter);
    return cursor.charFormat();
}

}  // namespace

class JobRunnerTest : public QObject {
    Q_OBJECT

private slots:
    void successDeliversTheResultOnTheInterfaceThread() {
        Lines log;
        JobRunner runner(log.function());
        QCOMPARE(runner.statusText(), QStringLiteral("Idle"));
        QSignalSpy running(&runner, &JobRunner::runningChanged);
        QSignalSpy finished(&runner, &JobRunner::finished);
        int result = 0;
        QThread* deliveredOn = nullptr;
        std::thread::id jobThread;
        QVERIFY(runner.start(
            QStringLiteral("Summing"),
            [&jobThread](const Reporter& reporter) {
                jobThread = std::this_thread::get_id();
                reporter.log("adding");
                reporter.progress(0.5, "half way");
                return 42;
            },
            [&](int value) {
                result = value;
                deliveredOn = QThread::currentThread();
            }));
        QVERIFY(runner.running());
        QCOMPARE(runner.title(), QStringLiteral("Summing"));
        QCOMPARE(runner.statusText(), QStringLiteral("Summing …"));
        QTRY_COMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Succeeded);
        QCOMPARE(result, 42);
        QCOMPARE(deliveredOn, QThread::currentThread());
        QVERIFY(jobThread != std::this_thread::get_id());
        QVERIFY(!runner.running());
        QCOMPARE(runner.statusText(), QStringLiteral("Summing: done"));
        QCOMPARE(runner.fraction(), 1.0);
        QCOMPARE(running.size(), 2);
        QCOMPARE(log.get(), (QStringList{QStringLiteral("--- Summing ---"), QStringLiteral("adding")}));

        bool ran = false;  // void jobs
        QVERIFY(runner.start(QStringLiteral("Nothing"), [](const Reporter&) {}, [&ran] { ran = true; }));
        QTRY_VERIFY(ran);
    }

    void progressTextAndThrottle() {
        Lines log;
        JobRunner runner(log.function());
        QStringList texts;
        connect(&runner, &JobRunner::statusChanged, this,
                [&texts](double, const QString& text) { texts.append(text); });
        QSignalSpy finished(&runner, &JobRunner::finished);
        // 10 000 steps of 0.01 % over about 0.5 s: MDATK's throttle shows one per 0.5 % at most
        runner.start(
            QStringLiteral("Counting"),
            [](const Reporter& reporter) {
                for (int i = 1; i <= 10000; ++i) {
                    reporter.progress(i / 10000.0, "frame");
                    std::this_thread::sleep_for(50us);
                }
            },
            [] {});
        QTRY_COMPARE_WITH_TIMEOUT(finished.size(), 1, 20000);
        QVERIFY2(texts.size() <= 203, qPrintable(QString::number(texts.size())));
        QVERIFY2(texts.size() >= 12, qPrintable(QString::number(texts.size())));  // and it does show progress
        QCOMPARE(texts.front(), QStringLiteral("Counting …"));
        QCOMPARE(texts.back(), QStringLiteral("Counting: done"));
        qInfo("progress updates shown for 10000 reports: %d", static_cast<int>(texts.size()) - 2);

        texts.clear();
        std::atomic<bool> release{false};
        runner.start(
            QStringLiteral("Half"),
            [&release](const Reporter& reporter) {
                reporter.progress(0.5, "half way");
                reporter.progress(0.501);  // below the step: not shown
                while (!release) std::this_thread::sleep_for(1ms);
            },
            [] {});
        QTRY_COMPARE(runner.statusText(), QStringLiteral("  50%  half way"));
        QCOMPARE(runner.fraction(), 0.5);
        release = true;
        QTRY_COMPARE(finished.size(), 2);
    }

    void aSecondJobIsRefused() {
        Lines log;
        JobRunner runner(log.function());
        QSignalSpy errors(&runner, &JobRunner::errorReported);
        QSignalSpy finished(&runner, &JobRunner::finished);
        runner.start(
            QStringLiteral("First"),
            [](const Reporter& reporter) {
                for (;;) {
                    reporter.checkCancel();
                    std::this_thread::sleep_for(1ms);
                }
            },
            [] {});
        bool secondRan = false;
        QVERIFY(!runner.start(QStringLiteral("Second"), [&secondRan](const Reporter&) { secondRan = true; }, [] {}));
        QCOMPARE(errors.size(), 1);
        QCOMPARE(errors.at(0).at(0).toString(), QStringLiteral("Another job is still running."));
        QCOMPARE(errors.at(0).at(1).toString(), QString());
        QVERIFY(log.get().contains(QStringLiteral("ERROR: Another job is still running.")));
        QCOMPARE(runner.title(), QStringLiteral("First"));
        runner.cancel();
        QTRY_COMPARE(finished.size(), 1);
        QVERIFY(!secondRan);
    }

    void cancelBeforeAJob() {
        Lines log;
        JobRunner runner(log.function());
        QVERIFY(!runner.cancel());
        QCOMPARE(runner.statusText(), QStringLiteral("Idle"));
        QVERIFY(log.get().isEmpty());

        // cancelled at once: whether the worker had begun or not, the job ends "Cancelled" unrun
        QSignalSpy finished(&runner, &JobRunner::finished);
        QSignalSpy errors(&runner, &JobRunner::errorReported);
        std::atomic<int> steps{0};
        bool succeeded = false;
        runner.start(
            QStringLiteral("Early"),
            [&steps](const Reporter& reporter) {
                for (;;) {
                    reporter.checkCancel();
                    ++steps;
                    std::this_thread::sleep_for(1ms);
                }
            },
            [&succeeded] { succeeded = true; });
        QVERIFY(runner.cancel());
        QCOMPARE(runner.statusText(), QStringLiteral("Cancelling …"));
        QTRY_COMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Cancelled);
        QVERIFY(steps <= 1);
        QVERIFY(!succeeded);
        QCOMPARE(errors.size(), 0);
        QCOMPARE(runner.statusText(), QStringLiteral("Cancelled"));
        QCOMPARE(log.get().back(), QStringLiteral("Cancelled."));
    }

    void cancelDuringAJob() {
        Lines log;
        JobRunner runner(log.function());
        QSignalSpy finished(&runner, &JobRunner::finished);
        QSignalSpy errors(&runner, &JobRunner::errorReported);
        bool succeeded = false;
        runner.start(
            QStringLiteral("Long"),
            [](const Reporter& reporter) {
                for (int i = 0;; ++i) {
                    reporter.checkCancel();
                    reporter.progress(std::min(0.99, i / 1000.0), "working");
                    std::this_thread::sleep_for(1ms);
                }
            },
            [&succeeded] { succeeded = true; });
        QTRY_VERIFY(runner.fraction() >= 0.02);
        QVERIFY(runner.cancel());
        QTRY_COMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Cancelled);
        QVERIFY(!succeeded);
        QCOMPARE(errors.size(), 0);
        QCOMPARE(runner.statusText(), QStringLiteral("Cancelled"));
        QCOMPARE(log.get().back(), QStringLiteral("Cancelled."));
        // the runner takes a new job
        bool again = false;
        QVERIFY(runner.start(QStringLiteral("Again"), [](const Reporter&) { return 1; }, [&again](int) { again = true; }));
        QTRY_VERIFY(again);
    }

    void cancelAfterAJob() {
        Lines log;
        JobRunner runner(log.function());
        QSignalSpy finished(&runner, &JobRunner::finished);
        runner.start(QStringLiteral("Quick"), [](const Reporter&) { return 7; }, [](int) {});
        QTRY_COMPARE(finished.size(), 1);
        QVERIFY(!runner.cancel());
        QCOMPARE(runner.statusText(), QStringLiteral("Quick: done"));
        QVERIFY(!log.get().contains(QStringLiteral("Cancelled.")));
        QCOMPARE(finished.size(), 1);
    }

    void failuresReachTheInterface() {
        Lines log;
        JobRunner runner(log.function());
        QSignalSpy errors(&runner, &JobRunner::errorReported);
        QSignalSpy finished(&runner, &JobRunner::finished);
        std::vector<std::string> order;
        JobFailure seen;
        bool wasValueError = false;
        connect(&runner, &JobRunner::errorReported, this, [&order] { order.push_back("dialog"); });
        runner.start(
            QStringLiteral("Analysing"),
            [](const Reporter&) -> int { throw mdbridge::core::ValueError("bad selection 'resid x'"); },
            [](int) { QFAIL("no result expected"); },
            [&](const JobFailure& failure) {
                order.push_back("handler");
                seen = failure;
                try {
                    std::rethrow_exception(failure.error);
                } catch (const mdbridge::core::ValueError&) {
                    wasValueError = true;
                } catch (...) {
                }
            });
        QTRY_COMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Failed);
        QCOMPARE(seen.title, QStringLiteral("Analysing"));
        QCOMPARE(seen.message, QStringLiteral("bad selection 'resid x'"));
        QVERIFY(wasValueError);
        QCOMPARE(order, (std::vector<std::string>{"handler", "dialog"}));
        QCOMPARE(errors.at(0).at(0).toString(), QStringLiteral("Analysing failed"));
        QCOMPARE(errors.at(0).at(1).toString(), QStringLiteral("bad selection 'resid x'"));
        QCOMPARE(runner.statusText(), QStringLiteral("Failed"));
        QCOMPARE(log.get().back(), QStringLiteral("ERROR: Analysing failed: bad selection 'resid x'"));

        // an exception that is no std::exception, and a result handler that throws
        runner.start(QStringLiteral("Odd"), [](const Reporter&) { throw 3; }, [] {});
        QTRY_COMPARE(finished.size(), 2);
        QCOMPARE(errors.at(1).at(1).toString(), QStringLiteral("unknown error"));
        runner.start(QStringLiteral("Showing"), [](const Reporter&) { return 1; },
                     [](int) { throw std::runtime_error("cannot show the result"); });
        QTRY_COMPARE(finished.size(), 3);
        QCOMPARE(finished.at(2).at(0).value<JobOutcome>(), JobOutcome::Failed);
        QCOMPARE(errors.at(2).at(0).toString(), QStringLiteral("Showing failed"));
        QCOMPARE(errors.at(2).at(1).toString(), QStringLiteral("cannot show the result"));
    }

    void busyAndStoppableWidgets() {
        Lines log;
        JobRunner runner(log.function());
        QPushButton render, cancel;
        QAction stopScripts(nullptr);
        runner.addBusy(&render);
        runner.addStoppable(&cancel);
        runner.addStoppable(&stopScripts);
        QVERIFY(render.isEnabled());
        QVERIFY(!cancel.isEnabled());
        QVERIFY(!stopScripts.isEnabled());
        std::atomic<bool> release{false};
        QSignalSpy finished(&runner, &JobRunner::finished);
        runner.start(QStringLiteral("Wait"), [&release](const Reporter&) { while (!release) std::this_thread::sleep_for(1ms); }, [] {});
        QVERIFY(!render.isEnabled());
        QVERIFY(cancel.isEnabled());
        QVERIFY(stopScripts.isEnabled());
        release = true;
        QTRY_COMPARE(finished.size(), 1);
        QVERIFY(render.isEnabled());
        QVERIFY(!cancel.isEnabled());
        QVERIFY(!stopScripts.isEnabled());
    }

    // AP-14: closing during a job asks; on "close" the job is cancelled and its own clean-up is
    // waited for before the window may close
    void closeWithoutAJobDoesNotAsk() {
        Lines log;
        JobRunner runner(log.function());
        int asked = 0;
        runner.setAsker([&asked](QWidget*, CloseQuestion, const QString&) {
            ++asked;
            return CloseAnswer::Cancel;
        });
        QVERIFY(runner.confirmClose(nullptr));
        QCOMPARE(asked, 0);
    }

    void closeQuestionAnsweredCancelKeepsTheJob() {
        Lines log;
        JobRunner runner(log.function());
        auto stop = std::make_shared<std::atomic<bool>>(false);
        QSignalSpy finished(&runner, &JobRunner::finished);
        runner.start(QStringLiteral("Rendering"), [stop](const Reporter& reporter) {
            while (!*stop) {
                reporter.checkCancel();
                std::this_thread::sleep_for(1ms);
            }
        }, [] {});
        QString askedAbout;
        runner.setAsker([&](QWidget*, CloseQuestion question, const QString& title) {
            askedAbout = title;
            if (question != CloseQuestion::JobRunning) QTest::qFail("unexpected question", __FILE__, __LINE__);
            return CloseAnswer::Cancel;
        });
        QVERIFY(!runner.confirmClose(nullptr));
        QCOMPARE(askedAbout, QStringLiteral("Rendering"));
        QVERIFY(runner.running());
        QCOMPARE(runner.statusText(), QStringLiteral("Rendering …"));  // not cancelled
        *stop = true;
        QTRY_COMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Succeeded);
    }

    void closeWaitsForTheJobsCleanUp() {
        QTemporaryDir dir;
        const QString temporary = dir.filePath(QStringLiteral("video.tmp.mp4"));
        Lines log;
        JobRunner runner(log.function());
        auto cleaned = std::make_shared<std::atomic<bool>>(false);
        runner.start(QStringLiteral("Rendering"), [temporary, cleaned](const Reporter& reporter) {
            QFile file(temporary);
            if (!file.open(QIODevice::WriteOnly)) throw std::runtime_error("cannot write");
            try {
                for (;;) {
                    reporter.checkCancel();
                    file.write("frame");
                    std::this_thread::sleep_for(1ms);
                }
            } catch (const mdbridge::core::Cancelled&) {
                std::this_thread::sleep_for(300ms);  // the job's clean-up takes a while
                file.close();
                file.remove();
                reporter.log("removed the partial video");
                *cleaned = true;
                throw;
            }
        }, [] {});
        QTRY_VERIFY(QFile::exists(temporary));
        std::vector<CloseQuestion> questions;
        runner.setAsker([&questions](QWidget*, CloseQuestion question, const QString&) {
            questions.push_back(question);
            return CloseAnswer::Close;
        });
        QElapsedTimer clock;
        clock.start();
        QVERIFY(runner.confirmClose(nullptr));
        QVERIFY(clock.elapsed() >= 250);
        QVERIFY(*cleaned);
        QVERIFY(!QFile::exists(temporary));
        QVERIFY(!runner.running());
        QCOMPARE(questions, std::vector<CloseQuestion>{CloseQuestion::JobRunning});
        QCOMPARE(runner.statusText(), QStringLiteral("Cancelled"));
        QVERIFY(log.get().contains(QStringLiteral("removed the partial video")));
    }

    void closeAsksAgainWhenTheJobDoesNotStop() {
        auto log = std::make_unique<Lines>();
        auto runner = std::make_unique<JobRunner>(log->function());
        runner->setStopWait(150ms);
        auto release = std::make_shared<std::atomic<bool>>(false);
        auto ended = std::make_shared<std::atomic<bool>>(false);
        auto began = std::make_shared<std::atomic<bool>>(false);
        runner->start(QStringLiteral("Stubborn"), [release, ended, began](const Reporter&) {
            *began = true;
            while (!*release) std::this_thread::sleep_for(1ms);  // never looks at the cancel flag
            *ended = true;
        }, [] {});
        QTRY_VERIFY(*began);
        std::vector<CloseQuestion> questions;
        std::vector<CloseAnswer> answers{CloseAnswer::Close, CloseAnswer::KeepWaiting, CloseAnswer::Cancel};
        runner->setAsker([&](QWidget*, CloseQuestion question, const QString&) {
            questions.push_back(question);
            const CloseAnswer answer = answers.front();
            answers.erase(answers.begin());
            return answer;
        });
        QVERIFY(!runner->confirmClose(nullptr));  // close, keep waiting, then stay open
        QCOMPARE(questions, (std::vector<CloseQuestion>{CloseQuestion::JobRunning, CloseQuestion::JobNotStopped,
                                                        CloseQuestion::JobNotStopped}));
        QVERIFY(runner->running());

        // a close request while the question is shown is refused (no second question)
        bool nested = true;
        answers = {CloseAnswer::Close, CloseAnswer::Close};
        questions.clear();
        runner->setAsker([&](QWidget*, CloseQuestion question, const QString&) {
            questions.push_back(question);
            if (questions.size() == 1) nested = runner->confirmClose(nullptr);
            const CloseAnswer answer = answers.front();
            answers.erase(answers.begin());
            return answer;
        });
        QVERIFY(runner->confirmClose(nullptr));  // close anyway: abandoned
        QVERIFY(!nested);
        QCOMPARE(questions.size(), std::size_t{2});
        QVERIFY(log->get().back().startsWith(QStringLiteral("WARNING: closed while Stubborn")));
        runner.reset();  // does not wait for the abandoned job
        QVERIFY(!*ended);
        *release = true;
        QTRY_VERIFY(*ended);
    }

    void destroyingTheRunnerWaitsForTheJob() {
        Lines log;
        auto runner = std::make_unique<JobRunner>(log.function());
        auto cleaned = std::make_shared<std::atomic<bool>>(false);
        auto began = std::make_shared<std::atomic<bool>>(false);
        runner->start(QStringLiteral("Writing"), [cleaned, began](const Reporter& reporter) {
            *began = true;
            try {
                for (;;) {
                    reporter.checkCancel();
                    std::this_thread::sleep_for(1ms);
                }
            } catch (const mdbridge::core::Cancelled&) {
                std::this_thread::sleep_for(200ms);
                *cleaned = true;
                throw;
            }
        }, [] {});
        QTRY_VERIFY(*began);
        runner.reset();
        QVERIFY(*cleaned);
    }
};

class LogViewTest : public QObject {
    Q_OBJECT

private slots:
    void keepsFiveToSixThousandLines() {
        LogView log;
        for (int chunk = 0; chunk < 7; ++chunk) {
            for (int i = 0; i < 1000; ++i) log.feed(QStringLiteral("line %1").arg(chunk * 1000 + i));
            log.flush();
            QVERIFY(log.lineCount() <= LogView::MaxLines);
        }
        QCOMPARE(log.lineCount(), LogView::KeepLines);  // 7000 > 6000: down to 5000
        const QStringList lines = shownLines(log);
        QCOMPARE(lines.size(), LogView::KeepLines);
        QCOMPARE(lines.front(), QStringLiteral("line 2000"));
        QCOMPARE(lines.back(), QStringLiteral("line 6999"));
        QCOMPARE(log.textView()->document()->blockCount(), LogView::KeepLines);
    }

    void aFloodShowsTheNewestLinesAndANote() {
        LogView log;
        for (int i = 0; i < 5000; ++i) log.feed(QStringLiteral("line %1").arg(i));
        log.flush();
        const QStringList lines = shownLines(log);
        QVERIFY(lines.size() <= LogView::MaxBatch + 1);
        QCOMPARE(lines.back(), QStringLiteral("line 4999"));
        QVERIFY(lines.front().startsWith(QStringLiteral("… ")));
        QVERIFY(lines.front().endsWith(QStringLiteral(" line(s) not shown: they came faster than the log can show them")));
        const int skipped = lines.front().section(QLatin1Char(' '), 1, 1).toInt();
        QCOMPARE(skipped + lines.size() - 1, 5000);
        QCOMPARE(formatOfLine(log, 0).foreground().color(), QColor(QStringLiteral("#1a5fb4")));
    }

    void linesAndBatches() {
        LogView log;
        log.feed(QStringLiteral("one\ntwo\r\nthree\n"));
        log.feed(QString());
        QCOMPARE(log.lineCount(), 0);  // not yet: batched
        QTRY_COMPARE(log.lineCount(), 4);
        QCOMPARE(log.text(), QStringLiteral("one\ntwo\nthree\n\n"));
        QVERIFY(log.textView()->isReadOnly());
    }

    void styles() {
        LogView log;
        QPalette light = log.palette();
        light.setColor(QPalette::Base, Qt::white);
        log.setPalette(light);
        log.feed(QStringLiteral("ERROR: Rendering failed: no frames"));
        log.feed(QStringLiteral("WARNING: project.layout: ignoring 1 extra panel(s)"));
        log.feed(QStringLiteral("WAITING for step 3"));
        log.feed(QStringLiteral("cpptraj: warning on stderr"), LogView::Kind::Stderr);
        log.feed(QStringLiteral("==== step 03: done"), LogView::Kind::Status);
        log.feed(QStringLiteral("planner message"), LogView::Kind::Info);
        log.feed(QStringLiteral("plain output"));
        const QTextCharFormat error = formatOfLine(log, 0);
        QCOMPARE(error.foreground().color(), QColor(QStringLiteral("#c01c28")));
        QCOMPARE(error.fontWeight(), static_cast<int>(QFont::Bold));
        QCOMPARE(formatOfLine(log, 1).foreground().color(), QColor(QStringLiteral("#b56a00")));
        QCOMPARE(formatOfLine(log, 1).fontWeight(), static_cast<int>(QFont::Bold));
        QCOMPARE(formatOfLine(log, 2).foreground().color(), QColor(QStringLiteral("#b56a00")));
        QCOMPARE(formatOfLine(log, 3).foreground().color(), QColor(QStringLiteral("#a51d2d")));
        QCOMPARE(formatOfLine(log, 4).fontWeight(), static_cast<int>(QFont::Bold));
        QCOMPARE(formatOfLine(log, 5).foreground().color(), QColor(QStringLiteral("#1a5fb4")));
        QVERIFY(!formatOfLine(log, 6).hasProperty(QTextFormat::ForegroundBrush));
        QCOMPARE(log.lineCount(), 7);

        QPalette dark = log.palette();  // lighter shades on a dark base
        dark.setColor(QPalette::Base, QColor(30, 30, 30));
        log.setPalette(dark);
        log.feed(QStringLiteral("ERROR: again"));
        QCOMPARE(formatOfLine(log, 7).foreground().color(), QColor(QStringLiteral("#ff7b63")));
    }

    void clearCopySave() {
        LogView log;
        log.feed(QStringLiteral("first\nsecond"));
        QCOMPARE(log.text(), QStringLiteral("first\nsecond\n"));
        log.copy();
        QCOMPARE(QApplication::clipboard()->text(), QStringLiteral("first\nsecond\n"));

        QTemporaryDir dir;
        const QString path = dir.filePath(QStringLiteral("log.txt"));
        QCOMPARE(log.saveAs(path), path);
        QFile file(path);
        QVERIFY(file.open(QIODevice::ReadOnly));
        QCOMPARE(QString::fromUtf8(file.readAll()), QStringLiteral("first\nsecond\n"));
        QCOMPARE(shownLines(log).back(), QStringLiteral("Log saved to %1").arg(path));

        const QString bad = dir.filePath(QStringLiteral("missing/folder/log.txt"));
        QCOMPARE(log.saveAs(bad), QString());
        QVERIFY(shownLines(log).back().startsWith(QStringLiteral("ERROR: could not write %1: ").arg(bad)));

        log.feed(QStringLiteral("queued"));  // pending lines go too
        log.clearAction()->trigger();
        QCOMPARE(log.lineCount(), 0);
        QCOMPARE(log.text(), QString());
        QTest::qWait(LogView::FlushMs * 3);
        QCOMPARE(log.lineCount(), 0);
        QCOMPARE(log.copyAction()->text(), QStringLiteral("Copy"));
        QCOMPARE(log.saveAction()->text(), QStringLiteral("Save…"));
        QCOMPARE(log.clearAction()->text(), QStringLiteral("Clear"));
    }

    void foldsWithCtrlL() {
        QWidget window;
        auto* log = new LogView(&window);
        window.resize(600, 300);
        log->resize(600, 300);
        window.show();
        window.activateWindow();
        QVERIFY(QTest::qWaitForWindowExposed(&window));
        QSignalSpy expanded(log, &LogView::expandedChanged);
        QCOMPARE(log->foldAction()->shortcut(), QKeySequence(Qt::CTRL | Qt::Key_L));
        QVERIFY(log->expanded());
        QVERIFY(log->textView()->isVisible());
        log->foldAction()->trigger();
        QVERIFY(!log->expanded());
        QVERIFY(!log->textView()->isVisible());
        QCOMPARE(expanded.size(), 1);
        QCOMPARE(expanded.at(0).at(0).toBool(), false);
        log->setExpanded(true);
        QVERIFY(log->textView()->isVisible());
        if (QTest::qWaitForWindowActive(&window)) {
            QTest::keyClick(&window, Qt::Key_L, Qt::ControlModifier);
            QVERIFY(!log->expanded());
        } else {
            qInfo("the platform activates no window: Ctrl+L not sent");
        }
    }

    void followsTheEnd() {
        LogView log;
        log.resize(400, 200);
        log.show();
        QVERIFY(QTest::qWaitForWindowExposed(&log));
        for (int i = 0; i < 500; ++i) log.feed(QStringLiteral("line %1").arg(i));
        log.flush();
        QScrollBar* bar = log.textView()->verticalScrollBar();
        QVERIFY(bar->maximum() > 0);
        QCOMPARE(bar->value(), bar->maximum());
        log.setFollowing(false);
        bar->setValue(0);
        log.feed(QStringLiteral("more"));
        log.flush();
        QCOMPARE(bar->value(), 0);
    }

    void aSinkOutlivesTheView() {
        std::function<void(const QString&)> sink;
        {
            LogView log;
            sink = log.sink(LogView::Kind::Stderr);
            sink(QStringLiteral("while it lives"));
            QCOMPARE(shownLines(log), QStringList{QStringLiteral("while it lives")});
            sink(QStringLiteral("pending at the end"));
        }
        sink(QStringLiteral("after"));  // nothing happens
        QTest::qWait(LogView::FlushMs * 2);
    }

    // 10^6 lines from a worker as fast as it can write: the loop never stalls for 0.2 s, the view
    // stays within its bounds and ends with the newest line
    void aMillionLinesFromAThread() {
        LogView log;
        log.resize(800, 400);
        log.show();
        QVERIFY(QTest::qWaitForWindowExposed(&log));
        std::atomic<bool> done{false};
        const auto sink = log.sink();
        QElapsedTimer clock;
        clock.start();
        std::thread writer([&done, sink] {
            for (int i = 0; i < 1000000; ++i) sink(QStringLiteral("line %1 of the pipeline output").arg(i));
            done = true;
        });
        int largest = 0;
        QTimer watch;
        connect(&watch, &QTimer::timeout, this, [&] { largest = std::max(largest, log.lineCount()); });
        watch.start(1);
        const double stall = longestStallSeconds([&done] { return !done.load(); });
        writer.join();
        const double feeding = static_cast<double>(clock.elapsed()) / 1000.0;
        QTRY_VERIFY(shownLines(log).back() == QStringLiteral("line 999999 of the pipeline output"));
        qInfo("10^6 lines from a thread: %.2f s, longest event-loop stall %.3f s, at most %d lines shown",
              feeding, stall, largest);
        QVERIFY2(stall < 0.2, qPrintable(QString::number(stall)));
        QVERIFY(largest <= LogView::MaxLines);
        QVERIFY(log.lineCount() <= LogView::MaxLines);
    }

    // 10^6 lines fed on the interface thread in chunks of 10 000 between events
    void aMillionLinesFromTheInterfaceThread() {
        LogView log;
        log.resize(800, 400);
        log.show();
        QVERIFY(QTest::qWaitForWindowExposed(&log));
        int next = 0;
        QTimer feeder;
        connect(&feeder, &QTimer::timeout, this, [&] {
            for (int i = 0; i < 10000 && next < 1000000; ++i, ++next) log.feed(QStringLiteral("line %1").arg(next));
        });
        feeder.start(0);
        QElapsedTimer clock;
        clock.start();
        const double stall = longestStallSeconds([&next] { return next < 1000000; });
        feeder.stop();
        qInfo("10^6 lines on the interface thread: %.2f s, longest event-loop stall %.3f s",
              static_cast<double>(clock.elapsed()) / 1000.0, stall);
        QVERIFY2(stall < 0.2, qPrintable(QString::number(stall)));
        QCOMPARE(shownLines(log).back(), QStringLiteral("line 999999"));
        QVERIFY(log.lineCount() <= LogView::MaxLines);
    }

    // a job flooding the log through its reporter: Cancel stays reachable and works at once
    void aFloodingJobCanBeCancelled() {
        LogView log;
        log.resize(800, 400);
        log.show();
        JobRunner runner(log.sink());
        QSignalSpy finished(&runner, &JobRunner::finished);
        std::atomic<long> written{0};
        runner.start(QStringLiteral("Chatty"), [&written](const Reporter& reporter) {
            for (long i = 0; i < 1000000; ++i) {
                reporter.checkCancel();
                reporter.rawLog("cpptraj line " + std::to_string(i));
                ++written;
            }
        }, [] {});
        QElapsedTimer clock;
        clock.start();
        QTimer::singleShot(150, &runner, [&runner] { runner.cancel(); });
        const double stall = longestStallSeconds([&finished] { return finished.isEmpty(); });
        qInfo("flooding job: cancelled after %ld lines, %.3f s; longest stall %.3f s", written.load(),
              static_cast<double>(clock.elapsed()) / 1000.0, stall);
        QCOMPARE(finished.size(), 1);
        QCOMPARE(finished.at(0).at(0).value<JobOutcome>(), JobOutcome::Cancelled);
        QVERIFY2(stall < 0.2, qPrintable(QString::number(stall)));
        QVERIFY(clock.elapsed() < 2000);
        QTRY_COMPARE(shownLines(log).back(), QStringLiteral("Cancelled."));
        QVERIFY(log.lineCount() <= LogView::MaxLines);
    }
};

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    int status = 0;
    {
        JobRunnerTest test;
        status |= QTest::qExec(&test, argc, argv);
    }
    {
        LogViewTest test;
        status |= QTest::qExec(&test, argc, argv);
    }
    return status;
}

#include "test_jobs.moc"
