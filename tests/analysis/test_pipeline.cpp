// The analysis cache and the analyze pipeline (src/analysis/cache.cpp, pipeline.cpp) against what
// MDATK's Pipeline does (data/analysis/pipeline_cases.json, tests/parity/make_pipeline_cases.py):
// the analyses a project needs, cache reuse and invalidation over sequences of runs, frame
// selections, external data, errors, a reused cleaned trajectory; the manifest's contents (keys,
// labels, units, time axis), the entries' descriptions and the dataset. Then MDBridge's own
// behaviour: the .npz format, cancellation, unreadable entries, MDATK's cache files left alone,
// cleaning without cleaned files, the frame after a run.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_test_support.hpp"
#include "mdbridge/analysis/cache.hpp"
#include "mdbridge/analysis/pipeline.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/hash.hpp"
#include "mdbridge/core/panels.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/version.hpp"

namespace mdbridge::test {
namespace {

namespace fs = std::filesystem;
using namespace analysis_cases;
using analysis::Pipeline;
using core::Table;

const Value& pipelineCases() { return caseFile("pipeline_cases.json"); }

void replaceAll(std::string& s, const std::string& from, const std::string& to) {
    if (from.empty()) {
        return;
    }
    for (std::size_t at = s.find(from); at != std::string::npos; at = s.find(from, at + to.size())) {
        s.replace(at, from.size(), to);
    }
}

std::string generic(const fs::path& path) {
    std::string s = core::pathToUtf8(path);
#ifdef _WIN32
    std::replace(s.begin(), s.end(), '\\', '/');
#endif
    return s;
}

// tests/data as MDBridge's file fingerprints spell it
const std::string& dataPrefix() {
    static const std::string prefix = [] {
        const std::string suffix = "/topology/demo_dry.parm7";
        std::string path = str(at(core::fileFingerprint(kData / "topology" / "demo_dry.parm7"), "path"));
#ifdef _WIN32
        std::replace(path.begin(), path.end(), '\\', '/');
#endif
        return path.substr(0, path.size() - suffix.size());
    }();
    return prefix;
}

// The generator's masking: {work} / {data} placeholders, cache keys of the masked fingerprints,
// mtime_ns 0; MDBridge's entry files named as MDATK's. On Windows every '\\' becomes '/' first
// (none of the compared texts holds one otherwise).
struct Masker {
    std::string work;
    std::map<std::string, std::string> keys;  // real -> fixed

    std::string text(std::string s) const {
#ifdef _WIN32
        // file fingerprints spell absolute paths with backslashes on Windows
        std::replace(s.begin(), s.end(), '\\', '/');
#endif
        replaceAll(s, work, "{work}");
        replaceAll(s, dataPrefix(), "{data}");
        for (const auto& [real, fixed] : keys) {
            replaceAll(s, real, fixed);
        }
        replaceAll(s, std::string(analysis::kEntryTag) + ".npz", ".npz");
        return s;
    }
    Value operator()(const Value& v) const {
        if (const Table* t = v.get<Table>()) {
            Table out;
            for (const auto& [k, x] : *t) {
                out.emplace_back(text(k), k == "mtime_ns" ? Value(0) : (*this)(x));
            }
            return out;
        }
        if (const core::Array* a = v.get<core::Array>()) {
            core::Array out;
            for (const Value& x : *a) {
                out.push_back((*this)(x));
            }
            return out;
        }
        if (const std::string* s = v.get<std::string>()) {
            return text(*s);
        }
        return v;
    }
    Value unmask(const Value& v) const {
        if (const Table* t = v.get<Table>()) {
            Table out;
            for (const auto& [k, x] : *t) {
                out.emplace_back(k, unmask(x));
            }
            return out;
        }
        if (const core::Array* a = v.get<core::Array>()) {
            core::Array out;
            for (const Value& x : *a) {
                out.push_back(unmask(x));
            }
            return out;
        }
        if (const std::string* s = v.get<std::string>()) {
            std::string r = *s;
            replaceAll(r, "{work}", work);
            replaceAll(r, "{data}", dataPrefix());
            return r;
        }
        return v;
    }
    std::vector<std::string> lines(const std::vector<std::string>& in) const {
        std::vector<std::string> out;
        for (const std::string& s : in) {
            out.push_back(text(s));
        }
        return out;
    }
};

// {type, message} as the generator records an exception (a KeyError by its argument)
Value errorRecord(const std::exception& e) {
    if (dynamic_cast<const std::out_of_range*>(&e)) {
        return Table{{"type", "KeyError"}, {"message", e.what()}};
    }
    const CaseError c = caseError(e);
    return Table{{"type", c.type}, {"message", c.message}};
}

std::string compact(const Value& v) { return core::toJson(sortedTables(untag(v)), -1); }

std::optional<std::vector<std::string>> namesOf(const Value& run) {
    const Value& names = at(run, "names");
    if (names.isNull()) {
        return std::nullopt;
    }
    return strings(names);
}

core::Array resultNames(const core::Dataset& d) {
    core::Array out;
    for (const core::AnalysisResult& r : d.results()) {
        out.push_back(r.name);
    }
    return out;
}

Value datasetRecord(const core::Dataset& d) {
    core::Array frames;
    for (const std::int64_t f : d.frames()) {
        frames.push_back(f);
    }
    Table labels;
    for (const std::string& key : d.seriesKeys()) {
        const core::Series s = d.series(key);
        labels.emplace_back(
            key, Table{{"label", s.label}, {"unit", s.unit}, {"target", s.target ? Value(*s.target) : Value()}});
    }
    for (const std::string& key : d.matrixKeys()) {
        const core::Matrix m = d.matrix(key);
        labels.emplace_back(key, Table{{"label", m.label}, {"kind", m.kind}, {"unit", m.unit},
                                       {"rows", static_cast<std::int64_t>(m.rowLabels.size())}});
    }
    for (const std::string& key : d.pairmapKeys()) {
        const core::PairMap p = d.pairmap(key);
        labels.emplace_back(key, Table{{"label", p.label},
                                       {"shape", core::Array{static_cast<std::int64_t>(p.shape().first),
                                                             static_cast<std::int64_t>(p.shape().second)}}});
    }
    return Table{{"frames", frames}, {"time_unit", d.timeUnit()}, {"results", resultNames(d)},
                 {"describe", d.describe()}, {"labels", labels}};
}

std::vector<std::string> cachedEntries(const std::string& dir) {
    std::vector<std::string> out;
    const std::string suffix = std::string(analysis::kEntryTag) + ".npz";
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(core::pathFromUtf8(dir), ec)) {
        const std::string name = core::pathToUtf8(entry.path().filename());
        if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            out.push_back(name.substr(0, name.size() - suffix.size()));
        }
    }
    return out;
}

// The keys of the masked fingerprint for every instance the run computed keys for.
void learnKeys(Masker& mask, Pipeline& pipe, const core::Project& project) {
    if (pipe.analysisKeys().empty()) {
        return;
    }
    const Value fp = mask(pipe.inputFingerprint(pipe.resolveInputs()));
    for (const auto& [instance, key] : pipe.analysisKeys()) {
        mask.keys[key] = analysis::analysisFor(instance).cacheKey(analysis::instanceProject(project, instance), fp);
    }
}

Value withoutKeys(Value v, std::initializer_list<std::string_view> drop) {
    Table* t = v.get<Table>();
    std::erase_if(*t, [&](const auto& e) { return std::find(drop.begin(), drop.end(), e.first) != drop.end(); });
    return v;
}

// ---- required_analyses -------------------------------------------------------------------------

TEST(Pipeline, RequiredAnalysesAsMdatk) {
    TempDir work;
    Masker mask{generic(work.path()), {}};
    for (const Value& c : list(at(pipelineCases(), "required"))) {
        const std::string name = str(at(c, "name"));
        SCOPED_TRACE(name);
        const core::Project project = core::projectFromPlain(mask.unmask(untag(at(c, "project"))));
        std::vector<std::string> logs;
        Pipeline pipe(project, core::Reporter([&](const std::string& m) { logs.push_back(m); }));
        const auto flag = [&](std::string_view key) -> std::optional<bool> {
            const Value& v = at(c, key);
            return v.isNull() ? std::nullopt : std::optional<bool>(*v.get<bool>());
        };
        try {
            const std::vector<std::string> names = pipe.requiredAnalyses(flag("video"), flag("charts"));
            ASSERT_TRUE(c.find("names")) << "MDATK failed: " << compact(at(c, "error"));
            EXPECT_EQ(strings(at(c, "names")), names);
        } catch (const std::exception& e) {
            ASSERT_TRUE(c.find("error")) << e.what();
            EXPECT_EQ(compact(at(c, "error")), compact(errorRecord(e)));
        }
        std::sort(logs.begin(), logs.end());
        EXPECT_EQ(strings(at(c, "log")), logs);
    }
}

TEST(Pipeline, PanelTypes) {
    std::vector<std::string> names;
    for (const core::PanelType& t : core::panelTypes()) {
        names.emplace_back(t.name);
    }
    EXPECT_EQ((std::vector<std::string>{"structure", "timeseries", "dual_axis", "heatmap", "contact_map", "sequence",
                                        "secondary"}),
              names);
}

// ---- analyze sequences ---------------------------------------------------------------------------

void prepareCleaned(const Value& seq, Masker& mask, const Value& firstRun) {
    for (const auto& [target, source] : *at(seq, "cleaned_files").get<Table>()) {
        const fs::path to = core::pathFromUtf8(str(mask.unmask(target)));
        fs::create_directories(to.parent_path());
        fs::copy_file(core::pathFromUtf8(str(mask.unmask(source))), to, fs::copy_options::overwrite_existing);
    }
    const core::Project project = core::projectFromPlain(mask.unmask(untag(at(firstRun, "project"))));
    const std::string key =
        analysis::cleaningKey(project.cleaning, project.inputs.topology, project.inputs.trajectories);
    writeFile(core::pathFromUtf8(str(mask.unmask(at(seq, "marker")))), "{\"key\": \"" + key + "\"}");
    // the key MDATK computes for masked fingerprints
    const std::string fixed = core::stableHash(
        mask(analysis::cleaningKeyInputs(project.cleaning, project.inputs.topology, project.inputs.trajectories)));
    EXPECT_EQ(str(at(seq, "cleaning_key")), fixed);
    mask.keys[key] = fixed;
}

void checkRun(const Value& run, Masker& mask, std::size_t index) {
    SCOPED_TRACE("run " + std::to_string(index) + " " + compact(at(run, "overrides")));
    const core::Project project = core::projectFromPlain(mask.unmask(untag(at(run, "project"))));
    std::vector<std::string> logs;
    std::vector<std::pair<double, std::string>> progress;
    Pipeline pipe(project, core::Reporter([&](const std::string& m) { logs.push_back(m); },
                                          [&](double f, const std::string& m) { progress.emplace_back(f, m); }));
    std::optional<core::Dataset> dataset;
    try {
        if (*at(run, "check").get<bool>()) {
            pipe.check();
        }
        dataset = pipe.analyze(namesOf(run), *at(run, "force").get<bool>());
    } catch (const std::exception& e) {
        try {
            learnKeys(mask, pipe, project);
        } catch (const std::exception&) {
        }
        ASSERT_TRUE(run.find("error")) << "MDBridge failed: " << e.what();
        EXPECT_EQ(compact(at(run, "error")), compact(mask(errorRecord(e))));
    }
    if (dataset) {
        ASSERT_FALSE(run.find("error")) << "MDATK failed: " << compact(at(run, "error"));
        learnKeys(mask, pipe, project);
    }
    EXPECT_EQ(strings(at(run, "log")), mask.lines(logs));
    const core::Array& expectedProgress = list(at(run, "progress"));
    ASSERT_EQ(expectedProgress.size(), progress.size());
    for (std::size_t n = 0; n < progress.size(); ++n) {
        EXPECT_EQ(hexFloat(list(expectedProgress[n])[0]), progress[n].first) << n;
        EXPECT_EQ(str(list(expectedProgress[n])[1]), mask.text(progress[n].second)) << n;
    }
    if (!dataset) {
        return;
    }
    // keys: the masked fingerprint's keys equal MDATK's
    core::Array keys;
    for (const auto& [instance, key] : pipe.analysisKeys()) {
        keys.push_back(core::Array{instance, mask.text(key)});
    }
    EXPECT_EQ(compact(at(run, "keys")), compact(keys));

    // the manifest
    const Value manifest = readJson(core::pathFromUtf8(pipe.cache().manifestPath()));
    std::vector<std::string> order;
    for (const auto& [k, v] : *manifest.get<Table>()) {
        order.push_back(k == "mdbridge_version" ? "mdatk_version" : k);
    }
    EXPECT_EQ(strings(at(run, "manifest_keys")), order);
    EXPECT_EQ(str(at(manifest, "project")), "MDBridge");
    EXPECT_EQ(str(at(manifest, "mdbridge_version")), kVersion);
    Value masked = mask(withoutKeys(manifest, {"project", "mdbridge_version", "updated"}));
    for (auto& [k, v] : *masked.get<Table>()) {
        if (k == "all_entries") {  // sorted by the real keys: by the masked ones, as the generator does
            core::Array& all = *v.get<core::Array>();
            const auto id = [](const Value& e) { return std::pair(str(at(e, "name")), str(at(e, "key"))); };
            std::stable_sort(all.begin(), all.end(), [&](const Value& x, const Value& y) { return id(x) < id(y); });
        }
    }
    EXPECT_EQ(compact(at(run, "manifest")), compact(masked));

    // the entries' descriptions
    for (const Value& entry : list(at(manifest, "entries"))) {
        const std::string name = str(at(entry, "name"));
        const Value info = readJson(core::pathFromUtf8(pipe.cache().entryBase(name, str(at(entry, "key"))) + ".json"));
        EXPECT_EQ(str(at(info, "project")), "MDBridge");
        EXPECT_EQ(str(at(info, "mdbridge_version")), kVersion);
        EXPECT_EQ(compact(at(at(run, "entries"), name)),
                  compact(mask(withoutKeys(info, {"project", "mdbridge_version", "created"}))))
            << name;
    }
    std::vector<std::string> cached;
    for (const std::string& e : cachedEntries(pipe.cache().cacheDir())) {
        cached.push_back(mask.text(e));
    }
    std::sort(cached.begin(), cached.end());
    EXPECT_EQ(strings(at(run, "cached")), cached);

    // the dataset
    const Value& expected = at(run, "dataset");
    std::vector<double> times;
    for (const Value& t : list(at(expected, "time_ps"))) {
        times.push_back(hexFloat(t));
    }
    EXPECT_EQ(times, dataset->timePs());
    EXPECT_EQ(compact(withoutKeys(expected, {"time_ps"})), compact(mask(datasetRecord(*dataset))));

    // and as load_cached_dataset rebuilds it
    const Value& reloaded = at(run, "reloaded");
    try {
        const core::Dataset again = Pipeline(project).loadCachedDataset();
        ASSERT_FALSE(reloaded.find("error")) << compact(reloaded);
        core::Array frames;
        for (const std::int64_t f : again.frames()) {
            frames.push_back(f);
        }
        EXPECT_EQ(compact(reloaded),
                  compact(Table{{"frames", frames}, {"describe", again.describe()}, {"results", resultNames(again)}}));
        // the stored results are the computed ones (the energy target is a display override)
        for (const core::AnalysisResult& r : again.results()) {
            const core::AnalysisResult* fresh = dataset->result(r.name);
            ASSERT_TRUE(fresh) << r.name;
            EXPECT_EQ(fresh->arrays, r.arrays) << r.name;
            EXPECT_EQ(core::toJson(fresh->params, -1), core::toJson(r.params, -1)) << r.name;
            EXPECT_EQ(core::toJson(fresh->meta, -1), core::toJson(r.meta, -1)) << r.name;
            EXPECT_EQ(fresh->key, r.key) << r.name;
        }
    } catch (const std::exception& e) {
        ASSERT_TRUE(reloaded.find("error")) << e.what();
        EXPECT_EQ(compact(at(reloaded, "error")), compact(mask(errorRecord(e))));
    }
}

class PipelineSequence : public ::testing::TestWithParam<std::string> {};

TEST_P(PipelineSequence, AsMdatk) {
    const Value* seq = nullptr;
    for (const Value& s : list(at(pipelineCases(), "sequences"))) {
        if (str(at(s, "name")) == GetParam()) {
            seq = &s;
        }
    }
    ASSERT_TRUE(seq);
    TempDir work;
    Masker mask{generic(work.path()), {}};
    for (const auto& [rel, text] : *at(*seq, "files").get<Table>()) {
        writeFile(work.path() / core::pathFromUtf8(rel), str(text));
    }
    const core::Array& runs = list(at(*seq, "runs"));
    if (seq->find("cleaned_files")) {
        prepareCleaned(*seq, mask, runs.front());
    }
    for (std::size_t n = 0; n < runs.size(); ++n) {
        checkRun(runs[n], mask, n);
        if (HasFatalFailure()) {
            return;
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Sequences, PipelineSequence,
                         ::testing::Values("reuse", "frames_external", "errors", "cleaned"));

// ---- MDBridge's own behaviour ------------------------------------------------------------------

// The demo project of the generator's base_project(), on its own cache folder.
core::Project demoProject(const TempDir& work) {
    const Value& first = list(at(list(at(pipelineCases(), "sequences")).front(), "runs")).front();
    Masker mask{generic(work.path()), {}};
    return core::projectFromPlain(mask.unmask(untag(at(first, "project"))));
}

TEST(Cache, NpzRoundTrip) {
    TempDir work;
    core::NamedArrays arrays;
    arrays.emplace_back("flags", core::NumArray(std::vector<bool>{true, false, true}));
    arrays.emplace_back("classes", core::NumArray(std::vector<std::int8_t>{-1, 0, 3, 127}, {2, 2}));
    arrays.emplace_back("map_i", core::NumArray(std::vector<std::int32_t>{1, -2, 3}));
    arrays.emplace_back("map_offsets", core::NumArray(std::vector<std::int64_t>{0, 1, 3, 3}));
    arrays.emplace_back("coords", core::NumArray(std::vector<float>{1.5f, -0.0f, 3.25f, 1e30f, 0.1f, 7.0f}, {1, 2, 3}));
    const double notANumber = std::numeric_limits<double>::quiet_NaN();
    arrays.emplace_back("rmsd_0", core::NumArray(std::vector<double>{0.1, notANumber, -2.5}));
    arrays.emplace_back("empty", core::NumArray(std::vector<double>{}, {0, 3}));
    arrays.emplace_back("resname", core::NumArray::strings({"G", "Ψψ", "", "LONGER"}, 4));
    arrays.emplace_back("chain_name", core::NumArray(std::vector<std::string>{"A", "B\xC3\xA4"}));
    std::vector<double> big(300000);
    for (std::size_t n = 0; n < big.size(); ++n) {
        big[n] = static_cast<double>(n % 977) * 0.5;
    }
    arrays.emplace_back("big", core::NumArray(big, {1000, 300}));
    const fs::path path = work / "a.npz";
    analysis::writeNpz(path, arrays);
    const core::NamedArrays back = analysis::readNpz(path);
    ASSERT_EQ(arrays.size(), back.size());
    for (std::size_t n = 0; n < arrays.size(); ++n) {
        EXPECT_EQ(arrays[n].first, back[n].first);
        EXPECT_EQ(arrays[n].second.dtypeName(), back[n].second.dtypeName()) << arrays[n].first;
        EXPECT_EQ(arrays[n].second.shape(), back[n].second.shape()) << arrays[n].first;
        if (arrays[n].first != "rmsd_0") {
            EXPECT_EQ(arrays[n].second, back[n].second) << arrays[n].first;
        }
    }
    const std::vector<double> rmsd = back[5].second.toDouble();
    EXPECT_EQ(0.1, rmsd[0]);
    EXPECT_TRUE(std::isnan(rmsd[1]));
    EXPECT_EQ(-2.5, rmsd[2]);
    EXPECT_EQ("<U4", std::string(back[7].second.dtypeName()));

    // damaged archives are errors
    std::string bytes = readFile(path);
    writeFile(work / "cut.npz", bytes.substr(0, bytes.size() / 2));
    EXPECT_THROW(analysis::readNpz(work / "cut.npz"), std::runtime_error);
    bytes[200] = static_cast<char>(bytes[200] ^ 0x55);
    writeFile(work / "flipped.npz", bytes);
    EXPECT_THROW(analysis::readNpz(work / "flipped.npz"), std::runtime_error);
}

TEST(Cache, ResultRoundTripAndUnreadableEntries) {
    TempDir work;
    analysis::CacheStore store(generic(work / "cache"));
    EXPECT_FALSE(store.has("rmsd", "0123456789abcdef"));
    EXPECT_FALSE(store.loadResult("rmsd", "0123456789abcdef"));
    core::AnalysisResult r;
    r.name = "contacts@A+B";
    r.key = "0123456789abcdef";
    r.params = Table{{"enabled", Value()}, {"cutoff", 4.5}};
    r.meta = Table{{"labels_a", core::Array{"G1", "C2"}},
                   {"label", "Pairwise contacts"},
                   {"x", std::numeric_limits<double>::infinity()}};
    r.addSeries("count", "count", "A–B contacts", "count");
    r.arrays.emplace_back("count", core::NumArray(std::vector<std::int64_t>{3, 4}));
    const std::vector<std::int64_t> frames{0, 5};
    const std::vector<double> times{0.0, 50.0};
    const std::string npz = store.saveResult(r, frames, times);
    EXPECT_EQ(store.entryBase(r.name, r.key) + ".npz", npz);
    EXPECT_TRUE(fs::exists(work / "cache" / "contacts@A+B-0123456789abcdef.mdbridge.npz"));
    EXPECT_TRUE(fs::exists(work / "cache" / "contacts@A+B-0123456789abcdef.mdbridge.json"));
    const auto loaded = store.loadResult(r.name, r.key);
    ASSERT_TRUE(loaded);
    EXPECT_EQ(frames, loaded->frames);
    EXPECT_EQ(times, loaded->timePs);
    EXPECT_EQ(r.arrays, loaded->result.arrays);
    EXPECT_EQ(core::toJson(r.meta, -1), core::toJson(loaded->result.meta, -1));
    EXPECT_EQ(core::toJson(r.params, -1), core::toJson(loaded->result.params, -1));

    // a damaged entry is ignored with a warning
    writeFile(core::pathFromUtf8(npz), "not a zip");
    LogCapture log;
    EXPECT_FALSE(store.loadResult(r.name, r.key));
    ASSERT_EQ(1u, log.lines().size());
    EXPECT_NE(std::string::npos,
              log.lines()[0].find("ignoring unreadable cache entry " + store.entryBase(r.name, r.key) + ": "));

    // no manifest
    try {
        store.loadDatasetFromManifest();
        ADD_FAILURE();
    } catch (const core::NamedError& e) {
        EXPECT_EQ("FileNotFoundError", e.pythonName());
        EXPECT_EQ("no MDBridge_cache.json in " + store.cacheDir(), std::string(e.what()));
    }
}

TEST(Pipeline, ReusesAndRecomputesDamagedEntriesButNotMdatkFiles) {
    TempDir work;
    const core::Project project = demoProject(work);
    std::vector<std::string> logs;
    const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); });
    Pipeline first(project, reporter);
    const core::Dataset d1 = first.analyze(std::vector<std::string>{"rmsd", "contacts"});
    const std::string rmsdKey = first.analysisKeys()[0].second;
    const std::string cache = first.cache().cacheDir();
    // MDATK's files for the same keys are never read
    EXPECT_FALSE(fs::exists(core::pathFromUtf8(cache) / "MDATK_cache.json"));
    writeFile(core::pathFromUtf8(cache) / ("rmsd-" + rmsdKey + ".npz"), "MDATK's");
    writeFile(core::pathFromUtf8(cache) / ("rmsd-" + rmsdKey + ".json"), "{}");
    // a damaged entry is computed again
    writeFile(core::pathFromUtf8(first.cache().entryBase("contacts", first.analysisKeys()[1].second) + ".npz"), "x");
    logs.clear();
    LogCapture log;
    Pipeline second(project, reporter);
    const core::Dataset d2 = second.analyze(std::vector<std::string>{"rmsd", "contacts"});
    EXPECT_EQ((std::vector<std::string>{"RMSD: loaded from cache (rmsd-" + rmsdKey + ")",
                                        "Preparing Pairwise contacts ...",
                                        "Analysing 24 frame(s): Pairwise contacts"}),
              logs);
    for (const core::AnalysisResult& r : d1.results()) {
        ASSERT_TRUE(d2.result(r.name));
        EXPECT_EQ(r.arrays, d2.result(r.name)->arrays) << r.name;
    }
    // the trajectory is back on its first frame
    EXPECT_EQ(0, second.system().trajectory->current().index);
}

TEST(Pipeline, CancellationStoresNothing) {
    TempDir work;
    const core::Project project = demoProject(work);
    auto cancel = std::make_shared<std::atomic<bool>>(false);
    std::vector<std::string> progress;
    const core::Reporter reporter({}, [&](double, const std::string& m) {
        progress.push_back(m);
        if (m == "Analysing frame 2/24") {
            cancel->store(true);
        }
    }, cancel);
    Pipeline pipe(project, reporter);
    EXPECT_THROW(pipe.analyze(std::vector<std::string>{"rmsd", "hbonds"}), core::Cancelled);
    EXPECT_EQ("Analysing frame 2/24", progress.back());
    EXPECT_TRUE(cachedEntries(pipe.cache().cacheDir()).empty());
    EXPECT_FALSE(fs::exists(core::pathFromUtf8(pipe.cache().manifestPath())));
    // a cancelled reporter stops a later run too; a fresh one completes it
    EXPECT_THROW(pipe.analyze(std::vector<std::string>{"rmsd"}), core::Cancelled);
    Pipeline again(project);
    const core::Dataset d = again.analyze(std::vector<std::string>{"rmsd", "hbonds"});
    EXPECT_EQ(24u, d.frameCount());
    EXPECT_EQ(2u, cachedEntries(again.cache().cacheDir()).size());
}

TEST(Pipeline, CleaningWithoutCleanedFiles) {
    TempDir work;
    core::Project project = demoProject(work);
    project.cleaning.enabled = true;
    Pipeline pipe(project);
    const std::string key =
        analysis::cleaningKey(project.cleaning, project.inputs.topology, project.inputs.trajectories);
    try {
        pipe.analyze(std::vector<std::string>{"rmsd"});
        ADD_FAILURE();
    } catch (const analysis::PipelineError& e) {
        const std::string what = e.what();
        EXPECT_EQ(0u, what.find("cleaning with cpptraj is not available in MDBridge yet: no cleaned trajectory for "
                                "these inputs and settings in "))
            << what;
        EXPECT_NE(std::string::npos, what.find("cleaned"));
        EXPECT_NE(std::string::npos, what.find(key + "; turn it off (--no-clean or cleaning.enabled = false) and give "
                                                     "the cleaned files as --top / --traj"))
            << what;
    }
    // a marker for other settings does not count
    const fs::path out = work / "cleaned";
    fs::create_directories(out);
    fs::copy_file(kData / "topology" / "demo_dry.parm7", out / "MDATK_cleaned.parm7");
    fs::copy_file(kData / "io" / "trajectories" / "demo_dry.nc", out / "MDATK_cleaned.nc");
    writeFile(out / "MDATK_clean.json", "{\"key\": \"0000000000000000\"}");
    project.cleaning.output_dir = generic(out);
    EXPECT_THROW(Pipeline(project).analyze(std::vector<std::string>{"rmsd"}), analysis::PipelineError);
    project.inputs.topology = generic(work / "missing.parm7");
    try {
        Pipeline(project).analyze(std::vector<std::string>{"rmsd"});
        ADD_FAILURE();
    } catch (const core::NamedError& e) {
        EXPECT_EQ("CpptrajError", e.pythonName());
        EXPECT_EQ("topology not found: " + project.inputs.topology, std::string(e.what()));
    }
}

TEST(Pipeline, ManifestKeepsEarlierEntries) {
    TempDir work;
    core::Project project = demoProject(work);
    Pipeline(project).analyze(std::vector<std::string>{"rmsd"});
    project.analyses.rmsd.reference_frame = 3;
    Pipeline second(project);
    second.analyze(std::vector<std::string>{"rmsd"});
    const Value manifest = readJson(core::pathFromUtf8(second.cache().manifestPath()));
    EXPECT_EQ(1u, list(at(manifest, "entries")).size());
    EXPECT_EQ(2u, list(at(manifest, "all_entries")).size());
    // an unreadable manifest is replaced (MDATK fails there)
    writeFile(core::pathFromUtf8(second.cache().manifestPath()), "{");
    Pipeline(project).analyze(std::vector<std::string>{"rmsd"});
    EXPECT_EQ(1u, list(at(readJson(core::pathFromUtf8(second.cache().manifestPath())), "all_entries")).size());
}

}  // namespace
}  // namespace mdbridge::test
