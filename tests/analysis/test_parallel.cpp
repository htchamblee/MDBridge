// The frame-parallel runner (analysis.hpp runAnalyses, plan W4.1): every analysis case of the
// parity data run with 1, 2, 4 and 8 threads gives byte-identical arrays, the same metadata, logs,
// progress calls, frame times and errors; ordered analyses see the sequential call order;
// errors, reader failures and cancellation behave as in the sequential run; the thread count.

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/core/json.hpp"

using namespace mdbridge;
using namespace mdbridge::test::analysis_cases;
using analysis::Analysis;
using analysis::AnalysisContext;
using analysis::AnalysisDescriptor;
using analysis::RunOptions;

namespace {

// The raw bytes of an array (strings: each value and a separator), with its dtype and shape.
std::string arrayBytes(const core::NumArray& a) {
    std::string out = std::string(a.dtypeName()) + "|";
    for (const std::size_t s : a.shape()) {
        out += std::to_string(s) + ",";
    }
    out += "|";
    const auto raw = [&out](const auto values) {
        out.append(reinterpret_cast<const char*>(values.data()), values.size_bytes());
    };
    switch (a.dtype()) {
    case core::NumArray::DType::Bool: raw(a.values<std::uint8_t>()); break;
    case core::NumArray::DType::Int8: raw(a.values<std::int8_t>()); break;
    case core::NumArray::DType::Int32: raw(a.values<std::int32_t>()); break;
    case core::NumArray::DType::Int64: raw(a.values<std::int64_t>()); break;
    case core::NumArray::DType::Float32: raw(a.values<float>()); break;
    case core::NumArray::DType::Float64: raw(a.values<double>()); break;
    case core::NumArray::DType::String:
        for (const std::string& s : a.values<std::string>()) {
            out += s;
            out += '\0';
        }
        break;
    }
    return out;
}

std::string doubleBytes(const std::vector<double>& values) {
    return {reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double)};
}

void expectSameRun(const CaseRun& expected, const CaseRun& actual) {
    EXPECT_EQ(expected.logs, actual.logs);
    ASSERT_EQ(expected.progress.size(), actual.progress.size());
    for (std::size_t i = 0; i < expected.progress.size(); ++i) {
        EXPECT_EQ(std::memcmp(&expected.progress[i].first, &actual.progress[i].first, sizeof(double)), 0) << i;
        EXPECT_EQ(expected.progress[i].second, actual.progress[i].second) << i;
    }
    ASSERT_EQ(expected.error.has_value(), actual.error.has_value());
    if (expected.error) {
        EXPECT_EQ(expected.error->type, actual.error->type);
        EXPECT_EQ(expected.error->message, actual.error->message);
        return;
    }
    EXPECT_EQ(doubleBytes(expected.run.timePs), doubleBytes(actual.run.timePs));
    ASSERT_EQ(expected.run.results.size(), actual.run.results.size());
    for (std::size_t r = 0; r < expected.run.results.size(); ++r) {
        const core::AnalysisResult& e = expected.run.results[r];
        const core::AnalysisResult& a = actual.run.results[r];
        SCOPED_TRACE(e.name);
        EXPECT_EQ(e.name, a.name);
        EXPECT_EQ(e.key, a.key);
        EXPECT_EQ(core::toJson(core::Value(e.params), -1), core::toJson(core::Value(a.params), -1));
        EXPECT_EQ(core::toJson(core::Value(e.meta), -1), core::toJson(core::Value(a.meta), -1));
        ASSERT_EQ(e.arrays.size(), a.arrays.size());
        for (std::size_t n = 0; n < e.arrays.size(); ++n) {
            EXPECT_EQ(e.arrays[n].first, a.arrays[n].first);
            EXPECT_TRUE(arrayBytes(e.arrays[n].second) == arrayBytes(a.arrays[n].second)) << e.arrays[n].first;
        }
    }
}

}  // namespace

// Every case of the analysis parity data, every analysis of the case, 1 thread against 2, 4, 8.
TEST(ParallelRunner, CasesIdenticalForAnyThreadCount) {
    std::size_t compared = 0;
    for (const std::string file :
         {"micro_cases.json", "demo_dry_cases.json", "demo_dry_variants.json", "systems_cases.json"}) {
        const Value& data = caseFile(file);
        for (const Value& c : list(at(data, "cases"))) {
            const std::string name = str(at(c, "name"));
            SCOPED_TRACE(file + " / " + name);
            const std::vector<std::string> instances = strings(at(c, "instances"));
            const CaseRun reference = runCase(data, c, instances, RunOptions{1});
            for (const std::size_t threads : std::vector<std::size_t>{2, 4, 8}) {
                SCOPED_TRACE("threads " + std::to_string(threads));
                expectSameRun(reference, runCase(data, c, instances, RunOptions{threads}));
            }
            ++compared;
        }
    }
    EXPECT_GE(compared, 100u);
}

// ---------------------------------------------------------------- runner semantics

namespace {

class TestDescriptor final : public AnalysisDescriptor {
public:
    explicit TestDescriptor(std::string label) : AnalysisDescriptor({"test", std::move(label), 1, false, 0, {}}) {}
    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { return nullptr; }
};

// Records its calls in order (ordered: process() on the calling thread, frame by frame).
class Ordered final : public Analysis {
public:
    Ordered(const AnalysisDescriptor& d, AnalysisContext& context, std::vector<std::string>& calls)
        : Analysis(d, context), m_calls(&calls) {}
    void prepare(const core::Frame& first) override {
        m_calls->push_back("prepare " + label() + " " + std::to_string(first.index));
    }
    void process(std::size_t k, const core::Frame& frame) override {
        m_calls->push_back("process " + label() + " " + std::to_string(k) + " " + std::to_string(frame.index));
    }
    core::AnalysisResult finalize() override { return newResult(); }

private:
    std::vector<std::string>* m_calls;
};

// Concurrent: per-frame slots (frame index and x), an error or a cancel at chosen frames.
class Concurrent final : public Analysis {
public:
    using Analysis::Analysis;
    std::optional<std::size_t> failAt, failAlsoAt, cancelAt;
    std::atomic<std::size_t> processed{0};
    std::int64_t preparedWith = -1;
    std::atomic<bool> processedBeforePrepare{false};

    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        m_index.assign(frames.size(), -1);
        m_x.assign(frames.size(), 0.0f);
    }
    void prepare(const core::Frame& first) override {
        if (processed.load() != 0) {
            processedBeforePrepare = true;
        }
        preparedWith = first.index;
    }
    bool concurrentProcess() const override { return true; }
    void process(std::size_t k, const core::Frame& frame) override {
        if (preparedWith < 0) {
            processedBeforePrepare = true;
        }
        ++processed;
        if (cancelAt && k == *cancelAt) {
            context().reporter().cancel();
        }
        if ((failAt && k == *failAt) || (failAlsoAt && k == *failAlsoAt)) {
            throw std::runtime_error("failed at frame " + std::to_string(k));
        }
        m_index[k] = frame.index;
        m_x[k] = frame.positions.at(0);
    }
    core::AnalysisResult finalize() override {
        core::AnalysisResult res = newResult();
        res.arrays.emplace_back("index", core::NumArray(m_index));
        res.arrays.emplace_back("x", core::NumArray(m_x));
        return res;
    }

private:
    std::vector<std::int64_t> m_index;
    std::vector<float> m_x;
};

// One atom at x = 0.5 * frame.
core::System lineSystem(std::size_t frames) {
    std::vector<std::vector<double>> coords;
    for (std::size_t f = 0; f < frames; ++f) {
        coords.push_back({0.5 * static_cast<double>(f), 0.0, 0.0});
    }
    return makeMicroSystem({{"C1", 1, "R", "C"}}, coords);
}

struct Outcome {
    std::vector<std::string> logs;
    std::vector<std::pair<double, std::string>> progress;
    std::vector<std::string> calls;
    std::optional<analysis::RunResult> run;
    std::string error;  // "<kind>: <what>"
    std::size_t processed = 0;
    std::int64_t preparedWith = -1;
    bool processedBeforePrepare = false;
    std::int64_t currentFrame = -1;
};

struct Faults {
    std::optional<std::size_t> failAt, failAlsoAt, cancelAt;
};

Outcome runTest(std::size_t systemFrames, const std::vector<std::int64_t>& frames, std::size_t threads,
                const Faults& setup = {}) {
    Outcome out;
    core::System system = lineSystem(systemFrames);
    const core::Reporter reporter([&](const std::string& m) { out.logs.push_back(m); },
                                  [&](double f, const std::string& m) { out.progress.emplace_back(f, m); });
    AnalysisContext context(system, core::Project(), "", reporter);
    const TestDescriptor first("First"), second("Second"), third("Third");
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(std::make_unique<Ordered>(first, context, out.calls));
    auto concurrent = std::make_unique<Concurrent>(second, context);
    concurrent->failAt = setup.failAt;
    concurrent->failAlsoAt = setup.failAlsoAt;
    concurrent->cancelAt = setup.cancelAt;
    Concurrent* probe = concurrent.get();
    analyses.push_back(std::move(concurrent));
    analyses.push_back(std::make_unique<Ordered>(third, context, out.calls));
    try {
        out.run = analysis::runAnalyses(context, analyses, frames, nullptr, RunOptions{threads});
    } catch (const core::Cancelled& e) {
        out.error = std::string("Cancelled: ") + e.what();
    } catch (const std::out_of_range& e) {
        out.error = std::string("out_of_range: ") + e.what();
    } catch (const std::runtime_error& e) {
        out.error = std::string("runtime_error: ") + e.what();
    }
    out.processed = probe->processed.load();
    out.preparedWith = probe->preparedWith;
    out.processedBeforePrepare = probe->processedBeforePrepare.load();
    out.currentFrame = system.currentFrame() ? system.currentFrame()->index : -1;
    return out;
}

std::vector<std::int64_t> range(std::int64_t begin, std::int64_t end, std::int64_t step = 1) {
    std::vector<std::int64_t> out;
    for (std::int64_t f = begin; f < end; f += step) {
        out.push_back(f);
    }
    return out;
}

void expectSameOutcome(const Outcome& expected, const Outcome& actual) {
    EXPECT_EQ(expected.logs, actual.logs);
    EXPECT_EQ(expected.progress, actual.progress);
    EXPECT_EQ(expected.calls, actual.calls);
    EXPECT_EQ(expected.error, actual.error);
    EXPECT_EQ(expected.preparedWith, actual.preparedWith);
    EXPECT_FALSE(actual.processedBeforePrepare);
    ASSERT_EQ(expected.run.has_value(), actual.run.has_value());
    if (expected.run) {
        EXPECT_EQ(doubleBytes(expected.run->timePs), doubleBytes(actual.run->timePs));
        ASSERT_EQ(expected.run->results.size(), actual.run->results.size());
        for (std::size_t r = 0; r < expected.run->results.size(); ++r) {
            EXPECT_EQ(expected.run->results[r].name, actual.run->results[r].name);
            ASSERT_EQ(expected.run->results[r].arrays.size(), actual.run->results[r].arrays.size());
            for (std::size_t n = 0; n < expected.run->results[r].arrays.size(); ++n) {
                EXPECT_TRUE(arrayBytes(expected.run->results[r].arrays[n].second) ==
                            arrayBytes(actual.run->results[r].arrays[n].second));
            }
        }
    }
}

}  // namespace

TEST(ParallelRunner, OrderedAndConcurrentAnalyses) {
    for (const auto& frames : {range(0, 1), range(0, 2), range(3, 450, 2), range(0, 1000)}) {
        SCOPED_TRACE(frames.size());
        const Outcome reference = runTest(1000, frames, 1);
        EXPECT_EQ(reference.processed, frames.size());
        EXPECT_EQ(reference.preparedWith, frames.front());
        EXPECT_EQ(reference.calls.front(), "prepare First " + std::to_string(frames.front()));
        EXPECT_EQ(reference.calls[1], "prepare Third " + std::to_string(frames.front()));
        EXPECT_EQ(reference.calls.size(), 2 + 2 * frames.size());
        EXPECT_EQ(reference.currentFrame, frames.back());
        for (const std::size_t threads : std::vector<std::size_t>{2, 3, 4, 8, 16}) {
            SCOPED_TRACE(threads);
            const Outcome parallel = runTest(1000, frames, threads);
            expectSameOutcome(reference, parallel);
            EXPECT_EQ(parallel.currentFrame, frames.back());  // the trajectory on the last analysed frame
        }
    }
}

// The error of the first failing frame, after the progress of the frames before it, as sequentially.
TEST(ParallelRunner, FirstErrorWins) {
    const std::vector<std::int64_t> frames = range(0, 300);
    for (const Faults& setup : {Faults{std::size_t{13}, std::size_t{29}, std::nullopt},
                               Faults{std::size_t{250}, std::size_t{0}, std::nullopt},
                               Faults{std::size_t{299}, std::nullopt, std::nullopt}}) {
        const Outcome reference = runTest(300, frames, 1, setup);
        ASSERT_EQ(reference.error.rfind("runtime_error: failed at frame ", 0), 0u) << reference.error;
        for (const std::size_t threads : std::vector<std::size_t>{2, 4, 8}) {
            SCOPED_TRACE(threads);
            expectSameOutcome(reference, runTest(300, frames, threads, setup));
        }
    }
}

// A frame the reader cannot read: the frames before it are finished first, then its error.
TEST(ParallelRunner, ReaderErrorAfterEarlierFrames) {
    std::vector<std::int64_t> frames = range(0, 40);
    frames.push_back(999);
    frames.push_back(1);
    const Outcome reference = runTest(50, frames, 1);
    ASSERT_EQ(reference.error.rfind("out_of_range: ", 0), 0u) << reference.error;
    EXPECT_EQ(reference.progress.size(), 40u);
    for (const std::size_t threads : std::vector<std::size_t>{2, 4, 8}) {
        SCOPED_TRACE(threads);
        expectSameOutcome(reference, runTest(50, frames, threads));
    }
    // a failing analysis frame before the unreadable one wins
    const Faults failing{std::size_t{7}, std::nullopt, std::nullopt};
    const Outcome failed = runTest(50, frames, 1, failing);
    ASSERT_EQ(failed.error, "runtime_error: failed at frame 7");
    expectSameOutcome(failed, runTest(50, frames, 4, failing));
}

// Cancelled from inside the run: Cancelled is thrown soon, progress stays monotone.
TEST(ParallelRunner, CancelStopsPromptly) {
    const std::vector<std::int64_t> frames = range(0, 5000);
    for (const std::size_t threads : std::vector<std::size_t>{1, 2, 4, 8}) {
        SCOPED_TRACE(threads);
        const Outcome out = runTest(5000, frames, threads, Faults{std::nullopt, std::nullopt, std::size_t{100}});
        EXPECT_EQ(out.error, "Cancelled: operation cancelled");
        EXPECT_LE(out.processed, 101 + 3 * threads);  // at most the frames already handed out
        for (std::size_t i = 1; i < out.progress.size(); ++i) {
            EXPECT_LT(out.progress[i - 1].first, out.progress[i].first);
        }
        ASSERT_FALSE(out.progress.empty());
        EXPECT_LE(out.progress.back().first, 101.0 / 5000.0);
    }
    // cancelled before the run: nothing processed, as sequentially
    for (const std::size_t threads : std::vector<std::size_t>{1, 4}) {
        core::System system = lineSystem(10);
        const core::Reporter reporter;
        reporter.cancel();
        AnalysisContext context(system, core::Project(), "", reporter);
        const TestDescriptor d("C");
        std::vector<std::unique_ptr<Analysis>> analyses;
        analyses.push_back(std::make_unique<Concurrent>(d, context));
        const std::vector<std::int64_t> some = range(0, 10);
        EXPECT_THROW((void)analysis::runAnalyses(context, analyses, some, nullptr, RunOptions{threads}),
                     core::Cancelled);
        EXPECT_EQ(static_cast<Concurrent&>(*analyses.front()).processed.load(), 0u);
    }
}

namespace {

// Sets an environment variable for the scope (unset when `value` is null), restoring it after.
class ScopedEnvironment {
public:
    ScopedEnvironment(const char* name, const char* value) : m_name(name) {
        if (const char* old = std::getenv(name)) {
            m_old = old;
        }
        set(value);
    }
    ~ScopedEnvironment() { set(m_old ? m_old->c_str() : nullptr); }
    ScopedEnvironment(const ScopedEnvironment&) = delete;
    ScopedEnvironment& operator=(const ScopedEnvironment&) = delete;

private:
    void set(const char* value) const {
#ifdef _WIN32
        _putenv_s(m_name.c_str(), value ? value : "");
#else
        if (value) {
            setenv(m_name.c_str(), value, 1);
        } else {
            unsetenv(m_name.c_str());
        }
#endif
    }
    std::string m_name;
    std::optional<std::string> m_old;
};

}  // namespace

TEST(ParallelRunner, ThreadCount) {
    using analysis::analysisThreads;
    EXPECT_EQ(analysisThreads(1), 1u);
    EXPECT_EQ(analysisThreads(6), 6u);
    EXPECT_EQ(analysisThreads(100000), analysis::kMaxThreads);
    {
        const ScopedEnvironment env(analysis::kThreadsEnvironment, "3");
        EXPECT_EQ(analysisThreads(), 3u);
        EXPECT_EQ(analysisThreads(2), 2u);  // the option wins
    }
    {
        const ScopedEnvironment env(analysis::kThreadsEnvironment, " 12 ");
        EXPECT_EQ(analysisThreads(), 12u);
    }
    {
        const ScopedEnvironment env(analysis::kThreadsEnvironment, "100000");
        EXPECT_EQ(analysisThreads(), analysis::kMaxThreads);
    }
    std::size_t fallback = 0;
    {
        const ScopedEnvironment env(analysis::kThreadsEnvironment, nullptr);
        fallback = analysisThreads();
        EXPECT_GE(fallback, 1u);
        EXPECT_LE(fallback, analysis::kMaxDefaultThreads);
    }
    for (const char* bad : {"0", "-2", "four", "2x", ""}) {
        const ScopedEnvironment env(analysis::kThreadsEnvironment, bad);
        EXPECT_EQ(analysisThreads(), fallback) << bad;
    }
}
