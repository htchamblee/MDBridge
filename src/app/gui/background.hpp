#pragma once
// Work off the interface's thread: jobs run on Qt's global thread pool and hand their results back
// to the thread of a receiver object, unless the receiver is gone by then (its Background was
// destroyed): the result is then dropped. Jobs must not throw (catch and return the problem).

#include <functional>
#include <memory>
#include <mutex>
#include <utility>

#include <QMetaObject>
#include <QObject>
#include <QThreadPool>

namespace mdbridge::app {

class Background {
public:
    explicit Background(QObject* receiver) : m_target(std::make_shared<Target>()) { m_target->receiver = receiver; }
    ~Background() {
        const std::lock_guard<std::mutex> lock(m_target->mutex);
        m_target->receiver = nullptr;
    }
    Background(const Background&) = delete;
    Background& operator=(const Background&) = delete;

    // job() on a pool thread, then done(result) on the receiver's thread.
    template <class Job, class Done>
    void run(Job job, Done done) const {
        using Result = decltype(job());
        auto shared = std::make_shared<std::pair<Job, Done>>(std::move(job), std::move(done));
        QThreadPool::globalInstance()->start(std::function<void()>([target = m_target, shared] {
            auto result = std::make_shared<Result>(shared->first());
            post(target, [shared, result] { shared->second(std::move(*result)); });
        }));
    }

    // A function that runs `fn` on the receiver's thread, callable from any thread for as long as
    // this Background lives (later calls do nothing).
    std::function<void(std::function<void()>)> poster() const {
        return [target = m_target](std::function<void()> fn) { post(target, std::move(fn)); };
    }

private:
    struct Target {
        std::mutex mutex;
        QObject* receiver = nullptr;
    };

    template <class Fn>
    static void post(const std::shared_ptr<Target>& target, Fn fn) {
        const std::lock_guard<std::mutex> lock(target->mutex);
        if (target->receiver) {
            // the receiver's posted events are discarded if it is deleted before they run
            QMetaObject::invokeMethod(target->receiver, std::move(fn), Qt::QueuedConnection);
        }
    }

    std::shared_ptr<Target> m_target;
};

}  // namespace mdbridge::app
