// The single-pass runner (mdatk/analysis/base.py run_analyses), frame-parallel: see analysis.hpp.

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <optional>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "mdbridge/analysis/analysis.hpp"

namespace mdbridge::analysis {

namespace {

constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

// A positive integer (surrounding blanks allowed), else 0.
std::size_t positiveInteger(std::string_view text) {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t')) {
        text.remove_prefix(1);
    }
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t')) {
        text.remove_suffix(1);
    }
    std::size_t value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc() || end != text.data() + text.size()) {
        return 0;
    }
    return value;
}

// One analysed frame on its way through the workers: a copy of the frame (the reader's frame
// changes on the next seek) and, after the workers, the first concurrent analysis that threw.
struct FrameSlot {
    std::size_t k = 0;
    core::Frame frame;
    bool done = false;  // guarded by the pool's mutex
    std::size_t failedAt = kNone;  // index into the analyses
    std::exception_ptr error;
};

// Worker threads running `work` on submitted slots; `work` must not throw. Destroying the pool
// drops the slots not started and waits for the running ones.
class WorkerPool {
public:
    WorkerPool(std::size_t threads, std::function<void(FrameSlot&)> work) : m_work(std::move(work)) {
        m_threads.reserve(threads);
        try {
            for (std::size_t t = 0; t < threads; ++t) {
                m_threads.emplace_back([this] { loop(); });
            }
        } catch (...) {
            shutdown();
            throw;
        }
    }
    ~WorkerPool() { shutdown(); }
    WorkerPool(const WorkerPool&) = delete;
    WorkerPool& operator=(const WorkerPool&) = delete;

    void submit(FrameSlot& slot) {
        {
            const std::lock_guard lock(m_mutex);
            slot.done = false;
            m_queue.push_back(&slot);
        }
        m_ready.notify_one();
    }

    // Until `slot` (submitted) is done; its results are then visible to the caller.
    void wait(const FrameSlot& slot) {
        std::unique_lock lock(m_mutex);
        m_finished.wait(lock, [&slot] { return slot.done; });
    }

private:
    void loop() {
        for (;;) {
            FrameSlot* slot = nullptr;
            {
                std::unique_lock lock(m_mutex);
                m_ready.wait(lock, [this] { return m_stop || !m_queue.empty(); });
                if (m_stop) {
                    return;
                }
                slot = m_queue.front();
                m_queue.pop_front();
            }
            m_work(*slot);
            {
                const std::lock_guard lock(m_mutex);
                slot->done = true;
            }
            m_finished.notify_all();
        }
    }

    void shutdown() {
        {
            const std::lock_guard lock(m_mutex);
            m_stop = true;
            m_queue.clear();
        }
        m_ready.notify_all();
        for (std::thread& t : m_threads) {
            t.join();
        }
        m_threads.clear();
    }

    std::function<void(FrameSlot&)> m_work;
    std::mutex m_mutex;
    std::condition_variable m_ready;
    std::condition_variable m_finished;
    std::deque<FrameSlot*> m_queue;
    bool m_stop = false;
    std::vector<std::thread> m_threads;
};

struct Progress {
    const core::Reporter& rep;
    std::size_t n;
    std::size_t every;

    void after(std::size_t k) const {
        if (k % every == 0 || k == n - 1) {
            rep.progress(static_cast<double>(k + 1) / static_cast<double>(n),
                         "Analysing frame " + std::to_string(k + 1) + "/" + std::to_string(n));
        }
    }
};

// Per-run checks of the analysed frames, in frame order on the calling thread (before the frame's
// non-concurrent process() calls, so logs are those of the sequential run):
// AN-19: with use_pbc, the first frame without a usable cell is reported once;
// AN-20: frames with non-finite coordinates are reported once per run of consecutive analysed
// frames ("frame k: ..." / "frames k1-k2: up to N atoms ...").
class FrameChecks {
public:
    explicit FrameChecks(const AnalysisContext& context) : m_context(context) {}

    void check(const core::Frame& frame) {
        if (m_context.usePbc() && !m_pbcWarned && !usableCell(frame.box)) {
            m_pbcWarned = true;
            m_context.warn("use_pbc is on but trajectory frame " + std::to_string(frame.index) +
                           " has no usable periodic cell; distances are not periodic");
        }
        std::size_t bad = 0;
        for (std::size_t i = 0; i + 2 < frame.positions.size(); i += 3) {
            if (!std::isfinite(frame.positions[i]) || !std::isfinite(frame.positions[i + 1]) ||
                !std::isfinite(frame.positions[i + 2])) {
                ++bad;
            }
        }
        if (bad == 0) {
            flush();
            return;
        }
        if (m_rangeFrames == 0) {
            m_first = frame.index;
        }
        m_last = frame.index;
        ++m_rangeFrames;
        m_maxBad = std::max(m_maxBad, bad);
    }

    void flush() {
        if (m_rangeFrames == 0) {
            return;
        }
        const std::string atoms = m_maxBad == 1 ? "1 atom has" : std::to_string(m_maxBad) + " atoms have";
        if (m_rangeFrames == 1) {
            m_context.warn("frame " + std::to_string(m_first) + ": " + atoms +
                           " non-finite coordinates; they are left out of distance searches");
        } else {
            m_context.warn("frames " + std::to_string(m_first) + "-" + std::to_string(m_last) + ": up to " + atoms +
                           " non-finite coordinates; they are left out of distance searches");
        }
        m_rangeFrames = 0;
        m_maxBad = 0;
    }

private:
    // three positive finite lengths and three finite angles in (0, 180)
    static bool usableCell(const std::optional<std::array<float, 6>>& box) {
        if (!box) {
            return false;
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const float length = (*box)[i];
            const float angle = (*box)[i + 3];
            if (!(std::isfinite(length) && length > 0.0f && std::isfinite(angle) && angle > 0.0f && angle < 180.0f)) {
                return false;
            }
        }
        return true;
    }

    const AnalysisContext& m_context;
    bool m_pbcWarned = false;
    std::int64_t m_first = 0;
    std::int64_t m_last = 0;
    std::size_t m_rangeFrames = 0;
    std::size_t m_maxBad = 0;
};

void prepareAll(std::span<Analysis* const> analyses, const core::Frame& first) {
    for (Analysis* a : analyses) {
        a->prepare(first);
    }
}

void runSequential(std::span<Analysis* const> analyses, std::span<const std::int64_t> frames,
                   core::Trajectory& trajectory, const Progress& progress, FrameChecks& checks, RunResult& out) {
    for (std::size_t k = 0; k < frames.size(); ++k) {
        progress.rep.checkCancel();
        const core::Frame& frame = trajectory.seek(frames[k]);
        out.timePs[k] = frame.time;
        if (k == 0) {
            prepareAll(analyses, frame);
        }
        checks.check(frame);
        for (Analysis* a : analyses) {
            a->process(k, frame);
        }
        progress.after(k);
    }
}

void runParallel(std::span<Analysis* const> analyses, std::span<const std::int64_t> frames,
                 core::Trajectory& trajectory, const Progress& progress, std::size_t threads, FrameChecks& checks,
                 RunResult& out) {
    const core::Reporter& rep = progress.rep;
    std::vector<bool> concurrent;
    std::vector<std::size_t> concurrentIndex;
    for (std::size_t i = 0; i < analyses.size(); ++i) {
        concurrent.push_back(analyses[i]->concurrentProcess());
        if (concurrent.back()) {
            concurrentIndex.push_back(i);
        }
    }
    const std::size_t n = frames.size();
    // frames read ahead of the oldest one not retired: enough to keep every worker busy
    std::vector<FrameSlot> slots(std::min(n, 2 * threads));
    std::vector<FrameSlot*> free;
    for (FrameSlot& s : slots) {
        free.push_back(&s);
    }
    std::deque<FrameSlot*> pending;  // submitted, in frame order
    // declared after the slots: destroyed (workers joined) before them
    WorkerPool pool(threads, [&](FrameSlot& slot) {
        slot.failedAt = kNone;
        slot.error = nullptr;
        if (rep.cancelled()) {
            return;  // the retirement throws Cancelled
        }
        for (const std::size_t i : concurrentIndex) {
            try {
                analyses[i]->process(slot.k, slot.frame);
            } catch (...) {
                slot.failedAt = i;
                slot.error = std::current_exception();
                return;
            }
        }
    });
    // The oldest pending frame, as the sequential run finishes it: the analyses in order (the
    // workers' error where the sequential run would have thrown it), then progress.
    const auto retire = [&] {
        FrameSlot& slot = *pending.front();
        pool.wait(slot);
        rep.checkCancel();
        checks.check(slot.frame);
        for (std::size_t i = 0; i < analyses.size(); ++i) {
            if (!concurrent[i]) {
                analyses[i]->process(slot.k, slot.frame);
            } else if (slot.failedAt == i) {
                std::rethrow_exception(slot.error);
            }
        }
        progress.after(slot.k);
        pending.pop_front();
        free.push_back(&slot);
    };
    for (std::size_t k = 0; k < n; ++k) {
        if (free.empty()) {
            retire();
        }
        FrameSlot& slot = *free.back();
        std::exception_ptr readError;
        try {
            rep.checkCancel();
            const core::Frame& frame = trajectory.seek(frames[k]);
            out.timePs[k] = frame.time;
            if (k == 0) {
                prepareAll(analyses, frame);
            }
            slot.k = k;
            slot.frame = frame;
        } catch (const core::Cancelled&) {
            throw;
        } catch (...) {
            readError = std::current_exception();
        }
        if (readError) {
            // the sequential run finishes the earlier frames first
            while (!pending.empty()) {
                retire();
            }
            std::rethrow_exception(readError);
        }
        free.pop_back();
        pending.push_back(&slot);
        pool.submit(slot);
    }
    while (!pending.empty()) {
        retire();
    }
}

}  // namespace

std::size_t analysisThreads(std::size_t requested) {
    if (requested > 0) {
        return std::min(requested, kMaxThreads);
    }
    if (const char* text = std::getenv(kThreadsEnvironment)) {
        if (const std::size_t value = positiveInteger(text); value > 0) {
            return std::min(value, kMaxThreads);
        }
    }
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware == 0 ? 1 : std::min<std::size_t>(hardware, kMaxDefaultThreads);
}

RunResult runAnalyses(AnalysisContext& context, std::span<Analysis* const> analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter, const RunOptions& options) {
    const core::Reporter& rep = reporter ? *reporter : context.reporter();
    core::Trajectory* trajectory = context.system().trajectory.get();
    if (!trajectory) {
        throw std::logic_error("runAnalyses: the system has no trajectory");
    }
    const std::size_t n = frames.size();
    std::string labels;  // ', '.join(a.label for a in analyses)
    for (std::size_t i = 0; i < analyses.size(); ++i) {
        rep.log("Preparing " + analyses[i]->label() + " ...");
        // SEL-17 / AN-13: every setup starts on the first analysed frame (an earlier setup, e.g.
        // rmsd's, may have moved the trajectory), so selections do not depend on the setup order
        if (n > 0 && trajectory->current().index != frames[0]) {
            trajectory->seek(frames[0]);
        }
        analyses[i]->setup(frames);
        labels += (i ? ", " : "") + analyses[i]->label();
    }
    RunResult out;
    out.timePs.resize(n);
    const Progress progress{rep, n, std::max<std::size_t>(1, n / 200)};
    rep.log("Analysing " + std::to_string(n) + " frame(s): " + labels);
    const std::size_t threads = std::min(analysisThreads(options.threads), n);
    const bool anyConcurrent =
        std::any_of(analyses.begin(), analyses.end(), [](const Analysis* a) { return a->concurrentProcess(); });
    FrameChecks checks(context);
    if (threads > 1 && anyConcurrent) {
        runParallel(analyses, frames, *trajectory, progress, threads, checks, out);
    } else {
        runSequential(analyses, frames, *trajectory, progress, checks, out);
    }
    checks.flush();
    out.results.reserve(analyses.size());
    for (Analysis* a : analyses) {
        out.results.push_back(a->finalize());
    }
    return out;
}

RunResult runAnalyses(AnalysisContext& context, const std::vector<std::unique_ptr<Analysis>>& analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter, const RunOptions& options) {
    std::vector<Analysis*> pointers;
    pointers.reserve(analyses.size());
    for (const auto& a : analyses) {
        pointers.push_back(a.get());
    }
    return runAnalyses(context, pointers, frames, reporter, options);
}

}  // namespace mdbridge::analysis
