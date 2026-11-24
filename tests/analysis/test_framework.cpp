// The analysis framework (mdatk/analysis/base.py, __init__.py, common.py): registry, instance ids,
// labels, expected keys, referenced groups, cache keys, results, the single-pass runner and the
// common helpers, against what MDATK computes for the cases of
// tests/parity/make_analysis_framework_cases.py (data/analysis, see data/README.md); plus the
// counterpart of MDATK's tests/test_v002.py::test_instance_keys_and_labels.

#include <algorithm>
#include <array>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/analysis/analyses.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/common.hpp"
#include "mdbridge/analysis/context.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::analysis;
using core::Array;
using core::Table;
using core::Value;

namespace {

const Value& cases() {
    static const Value data = test::readJson(test::kData / "analysis" / "framework_cases.json");
    return data;
}

const Value& at(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    if (!v) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *v;
}
const Array& list(const Value& value) { return *value.get<Array>(); }
const std::string& str(const Value& value) { return *value.get<std::string>(); }
std::int64_t integer(const Value& value) { return *value.get<std::int64_t>(); }
double hexFloat(const Value& value) { return std::strtod(str(value).c_str(), nullptr); }

std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& v : list(value)) {
        out.push_back(str(v));
    }
    return out;
}

template <class T>
std::vector<T> ints(const Value& value) {
    std::vector<T> out;
    for (const Value& v : list(value)) {
        out.push_back(static_cast<T>(integer(v)));
    }
    return out;
}

Value expectedOutputsValue(const ExpectedOutputs& outputs) {
    const auto array = [](const std::vector<std::string>& items) { return Value(Array(items.begin(), items.end())); };
    return Table{{"series", array(outputs.series)}, {"matrices", array(outputs.matrices)},
                 {"pairmaps", array(outputs.pairmaps)}};
}

// describe() of the generator for one analysis on one project
void expectDescription(const AnalysisDescriptor& d, const core::Project& project, const Value& expected,
                       const std::string& context) {
    SCOPED_TRACE(context);
    EXPECT_EQ(d.enabled(project), *at(expected, "enabled").get<bool>());
    EXPECT_EQ(core::pythonRepr(Value(d.paramsFor(project))), str(at(expected, "params")));
    EXPECT_EQ(d.selectionTexts(project), strings(at(expected, "selection_texts")));
    EXPECT_EQ(expectedOutputsValue(d.expectedOutputs(project)), at(expected, "expected_outputs"));
    EXPECT_EQ(d.groups(project), strings(at(expected, "groups")));
    for (const Value& pair : list(at(expected, "labels"))) {
        EXPECT_EQ(d.outputLabel(project, str(list(pair)[0])), str(list(pair)[1])) << str(list(pair)[0]);
    }
    const Array& inputs = list(at(cases(), "inputs"));
    const Array& keys = list(at(expected, "cache_keys"));
    ASSERT_EQ(inputs.size(), keys.size());
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        EXPECT_EQ(d.cacheKey(project, inputs[i]), str(keys[i])) << "inputs " << i;
    }
    EXPECT_EQ(Value(referencedGroups(d.selectionTexts(project), project)), at(expected, "referenced_groups"));
}

// {"type", "message"} of an error: KeyError is std::out_of_range, ValueError core::ValueError
template <class F>
void expectError(F&& call, const Value& expected, const std::string& context) {
    const std::string& type = str(at(expected, "type"));
    try {
        call();
        ADD_FAILURE() << context << ": expected " << type;
    } catch (const core::ValueError& e) {
        EXPECT_EQ(type, "ValueError") << context;
        EXPECT_EQ(e.what(), str(at(expected, "message"))) << context;
    } catch (const std::out_of_range& e) {
        EXPECT_EQ(type, "KeyError") << context;
        EXPECT_EQ(e.what(), str(at(expected, "message"))) << context;
    }
}

// A topology of `residueOfAtom.size()` atoms named A<i> in residues R<r> (resid r + 1), one segment.
std::shared_ptr<core::Topology> makeTopology(const std::vector<std::int32_t>& residueOfAtom) {
    auto t = std::make_shared<core::Topology>();
    std::int32_t residues = 0;
    for (std::size_t i = 0; i < residueOfAtom.size(); ++i) {
        t->names.push_back(t->strings.intern("A" + std::to_string(i)));
        t->atomResidue.push_back(residueOfAtom[i]);
        residues = std::max(residues, residueOfAtom[i] + 1);
    }
    for (std::int32_t r = 0; r < residues; ++r) {
        t->resids.push_back(r + 1);
        t->resnames.push_back(t->strings.intern("R" + std::to_string(r)));
        t->residueSegment.push_back(0);
    }
    t->segids.push_back(t->strings.intern("SYS"));
    t->finish();
    return t;
}

}  // namespace

// ---------------------------------------------------------------- registry and projects

TEST(AnalysisFramework, Registry) {
    const auto& registry = registeredAnalyses();
    const Array& expected = list(at(cases(), "registry"));
    ASSERT_EQ(registry.size(), expected.size());
    for (std::size_t i = 0; i < registry.size(); ++i) {
        const AnalysisDescriptor& d = *registry[i];
        const Value& e = expected[i];
        SCOPED_TRACE(d.name());
        EXPECT_EQ(d.name(), str(at(e, "name")));
        EXPECT_EQ(d.label(), str(at(e, "label")));
        EXPECT_EQ(d.version(), integer(at(e, "version")));
        EXPECT_EQ(d.internal(), *at(e, "internal").get<bool>());
        EXPECT_EQ(d.groupSlots(), integer(at(e, "group_slots")));
        OutputLabels labels;
        for (const Value& pair : list(at(e, "output_labels"))) {
            labels.emplace_back(str(list(pair)[0]), str(list(pair)[1]));
        }
        EXPECT_EQ(d.outputLabels(), labels);
        EXPECT_EQ(findAnalysis(d.name()), &d);
        EXPECT_EQ(registryIndex(d.name()), static_cast<int>(i));
    }
    EXPECT_EQ(findAnalysis("bogus"), nullptr);
    EXPECT_EQ(findAnalysis("contacts@A+B"), nullptr);
    EXPECT_EQ(registryIndex("bogus"), -1);
}

TEST(AnalysisFramework, Projects) {
    for (const Value& c : list(at(cases(), "projects"))) {
        const std::string name = str(at(c, "name"));
        SCOPED_TRACE(name);
        const core::Project project = core::projectFromPlain(at(c, "plain"));
        for (const AnalysisDescriptor* d : registeredAnalyses()) {
            const Value& e = at(at(c, "analyses"), d->name());
            expectDescription(*d, project, e, name + " / " + d->name());

            // new_result(key, **meta) of an analysis object
            core::System system{makeTopology({0}), nullptr};
            AnalysisContext context(system, project);
            struct Probe final : Analysis {
                using Analysis::Analysis;
                void process(std::size_t, const core::Frame&) override {}
                core::AnalysisResult finalize() override { return {}; }
            };
            const Array& results = list(at(e, "new_result"));
            const Probe plain(*d, context);
            const core::AnalysisResult a = plain.newResult();
            const Probe named(*d, context, d->name() + "@X");
            const core::AnalysisResult b =
                named.newResult("0123abcd", Table{{"x", 1}, {"label", "old"}, {"z", Array{1.5}}});
            for (const auto& [got, expected] : {std::pair{&a, &results[0]}, std::pair{&b, &results[1]}}) {
                EXPECT_EQ(got->name, str(at(*expected, "name")));
                EXPECT_EQ(got->key, str(at(*expected, "key")));
                EXPECT_EQ(core::pythonRepr(Value(got->params)), str(at(*expected, "params")));
                EXPECT_EQ(core::pythonRepr(Value(got->meta)), str(at(*expected, "meta")));
            }
        }
        EXPECT_EQ(expectedOutputsValue(expectedKeys(project)), at(c, "expected_keys"));
        EXPECT_EQ(expectedOutputsValue(expectedKeys(project, false)), at(c, "expected_keys_enabled"));
        for (const Value& pair : list(at(c, "referenced_groups"))) {
            EXPECT_EQ(Value(referencedGroups(strings(list(pair)[0]), project)), list(pair)[1])
                << core::toJson(list(pair)[0], -1);
        }
        for (const Value& entry : list(at(c, "instances"))) {
            const std::string instance = str(at(entry, "instance"));
            SCOPED_TRACE(instance);
            const Instance parsed = parseInstance(instance);
            EXPECT_EQ(parsed.base, str(list(at(entry, "parse"))[0]));
            EXPECT_EQ(parsed.groups, strings(list(at(entry, "parse"))[1]));
            EXPECT_EQ(instanceId(parsed.base, parsed.groups), str(at(entry, "id")));
            if (const Value* error = entry.find("error")) {
                expectError([&] { (void)instanceProject(project, instance); }, *error, instance);
                continue;
            }
            const core::Project q = instanceProject(project, instance);
            if (*at(entry, "unchanged").get<bool>()) {
                EXPECT_EQ(q, project);
            }
            expectDescription(analysisFor(instance), q, at(entry, "project"), instance);
        }
        for (const Value& pair : list(at(c, "key_labels"))) {
            EXPECT_EQ(keyLabel(str(list(pair)[0]), project), str(list(pair)[1])) << str(list(pair)[0]);
        }
    }
}

TEST(AnalysisFramework, InstanceIdsAndKeys) {
    for (const Value& c : list(at(cases(), "instance_id"))) {
        EXPECT_EQ(instanceId(str(list(c)[0]), strings(list(c)[1])), str(list(c)[2]));
    }
    for (const Value& c : list(at(cases(), "split_key"))) {
        const auto [instance, out] = splitKey(str(list(c)[0]));
        EXPECT_EQ(instance, str(list(list(c)[1])[0]));
        EXPECT_EQ(out, str(list(list(c)[1])[1]));
    }
}

// tests/test_v002.py::test_instance_keys_and_labels
TEST(AnalysisFramework, InstanceKeysAndLabels) {
    const core::Project p;
    EXPECT_EQ(parseInstance("contacts@Target+Aptamer"), (Instance{"contacts", {"Target", "Aptamer"}}));
    EXPECT_EQ(splitKey("hbonds@A+B/count"), (std::pair<std::string, std::string>{"hbonds@A+B", "count"}));
    const core::Project q = instanceProject(p, "contacts@Target+Aptamer");
    EXPECT_EQ(q.analyses.contacts.group_a, "Target");
    EXPECT_EQ(q.analyses.contacts.group_b, "Aptamer");
    EXPECT_EQ(p.analyses.contacts.group_a, "Aptamer");  // the original project is untouched
    EXPECT_EQ(keyLabel("contacts@Target+Aptamer/atom_pairs", p), "Target–Aptamer heavy-atom contacts");
    EXPECT_EQ(keyLabel("rmsd@Target/Backbone", p), "RMSD Backbone");
    EXPECT_THROW((void)instanceProject(p, "contacts@Nope+Target"), core::ValueError);
    EXPECT_THROW((void)instanceProject(p, "contacts@Target"), core::ValueError);
    EXPECT_THROW((void)instanceProject(p, "trace@Target"), core::ValueError);
    EXPECT_THROW((void)analysisFor("nope"), std::out_of_range);
    EXPECT_EQ(&analysisFor("rmsd@Target"), &rmsdAnalysis());
}

TEST(AnalysisFramework, RmsdOutputKeyAndStructureGroups) {
    EXPECT_EQ(rmsdOutputKey("Complex", 0), "Complex");
    EXPECT_EQ(rmsdOutputKey("  a  b ", 3), "a_b");
    EXPECT_EQ(rmsdOutputKey("x__ y", 0), "x___y");
    EXPECT_EQ(rmsdOutputKey("ÅÅ", 2), "target2");
    EXPECT_EQ(rmsdOutputKey("", 7), "target7");
    EXPECT_EQ(rmsdOutputKey("a.b+c-d", 0), "a.b+c-d");
    core::Project p;
    EXPECT_EQ(structureGroups(p), (std::vector<std::string>{"Aptamer", "Target"}));
    p.visual.structure_groups = {"Target", "Nope", "Target"};
    EXPECT_EQ(structureGroups(p), (std::vector<std::string>{"Target", "Target"}));
}

TEST(AnalysisFramework, ComputationsNotPortedYetThrow) {
    core::System system{makeTopology({0}), nullptr};
    AnalysisContext context(system, core::Project());
    for (const AnalysisDescriptor* d : registeredAnalyses()) {
        try {
            (void)d->create(context);
            // a ported analysis: nothing to check here (its own tests cover it)
        } catch (const core::NamedError& e) {
            EXPECT_EQ(e.pythonName(), "NotImplementedError");
            EXPECT_EQ(std::string(e.what()), d->name() + ": not implemented yet");
        }
    }
}

// ---------------------------------------------------------------- common.py

TEST(AnalysisFramework, ResidueRowsAndWithout) {
    const auto topology = makeTopology({0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 5, 5, 6, 7, 7, 7, 7});
    for (const Value& c : list(at(cases(), "residue_rows"))) {
        const core::AtomGroup atoms(topology, ints<std::int32_t>(at(c, "atoms")));
        const ResidueRows rows = residueRows(atoms);
        EXPECT_EQ(rows.residues, ints<std::int32_t>(at(c, "residues")));
        EXPECT_EQ(rows.rows, ints<std::int64_t>(at(c, "rows")));
    }
    for (const Value& c : list(at(cases(), "without"))) {
        const core::AtomGroup atoms(topology, ints<std::int32_t>(at(c, "atoms")));
        const core::AtomGroup other(topology, ints<std::int32_t>(at(c, "other")));
        EXPECT_EQ(without(atoms, other).indices(), ints<std::int32_t>(at(c, "result")));
    }
}

TEST(AnalysisFramework, PerRowCounts) {
    for (const Value& c : list(at(cases(), "per_row_counts"))) {
        const auto offsets = ints<std::int64_t>(at(c, "offsets"));
        const auto rows = ints<std::int32_t>(at(c, "rows"));
        std::vector<float> values;
        for (const Value& v : list(at(c, "values"))) {
            values.push_back(static_cast<float>(hexFloat(v)));
        }
        const auto nRows = static_cast<std::size_t>(integer(at(c, "n_rows")));
        const core::NumArray out = perRowCounts(offsets, rows, values, nRows);
        EXPECT_EQ(out.dtypeName(), "float32");
        EXPECT_EQ(out.shape(), ints<std::size_t>(at(c, "shape")));
        const auto got = out.values<float>();
        const Array& expected = list(at(c, "result"));
        ASSERT_EQ(got.size(), expected.size());
        for (std::size_t i = 0; i < got.size(); ++i) {
            EXPECT_EQ(static_cast<double>(got[i]), hexFloat(expected[i])) << i;
        }
    }
    const std::vector<std::int64_t> offsets{0, 1};
    const std::vector<std::int32_t> rows{2};
    const std::vector<float> values{1.0f};
    EXPECT_THROW((void)perRowCounts(offsets, rows, values, 2), std::out_of_range);
    EXPECT_THROW((void)perRowCounts({}, {}, {}, 2), std::invalid_argument);
    EXPECT_THROW((void)perRowCounts(offsets, {}, {}, 2), std::invalid_argument);
}

// ---------------------------------------------------------------- context and runner

namespace {

struct Recording {
    std::vector<std::string> calls;  // "setup L f0,f1", "process L k frame", "finalize L n"
};

class RecorderDescriptor final : public AnalysisDescriptor {
public:
    explicit RecorderDescriptor(std::string label)
        : AnalysisDescriptor({"recorder", std::move(label), 1, false, 0, {}}) {}
    std::unique_ptr<Analysis> create(AnalysisContext&, std::string) const override { return nullptr; }
};

class Recorder final : public Analysis {
public:
    Recorder(const AnalysisDescriptor& d, AnalysisContext& context, Recording& recording)
        : Analysis(d, context), m_recording(&recording) {}
    void setup(std::span<const std::int64_t> frames) override {
        Analysis::setup(frames);
        std::string text = "setup " + label() + " ";
        for (std::size_t i = 0; i < frames.size(); ++i) {
            text += (i ? "," : "") + std::to_string(frames[i]);
        }
        m_recording->calls.push_back(text);
    }
    void process(std::size_t k, const core::Frame& frame) override {
        m_recording->calls.push_back("process " + label() + " " + std::to_string(k) + " " +
                                     std::to_string(frame.index));
    }
    core::AnalysisResult finalize() override {
        m_recording->calls.push_back("finalize " + label() + " " + std::to_string(frameCount()));
        core::AnalysisResult result;
        result.name = label();
        return result;
    }

private:
    Recording* m_recording;
};

core::System memorySystem(std::size_t frames, double dt) {
    std::vector<core::Frame> list(frames);
    for (core::Frame& f : list) {
        f.positions = {0.0f, 0.0f, 0.0f};
    }
    return core::System{makeTopology({0}), std::make_unique<core::MemoryTrajectory>(std::move(list), dt)};
}

}  // namespace

TEST(AnalysisFramework, RunAnalyses) {
    for (const Value& c : list(at(cases(), "run"))) {
        const auto frames = ints<std::int64_t>(at(c, "frames"));
        const std::vector<std::string> labels = strings(at(c, "labels"));
        SCOPED_TRACE(core::toJson(at(c, "labels"), -1) + " " + std::to_string(frames.size()));
        core::System system = memorySystem(static_cast<std::size_t>(integer(at(c, "n_total"))), hexFloat(at(c, "dt")));
        std::vector<std::string> logs;
        std::vector<std::pair<double, std::string>> progress;
        const core::Reporter reporter([&](const std::string& m) { logs.push_back(m); },
                                      [&](double f, const std::string& m) { progress.emplace_back(f, m); });
        AnalysisContext context(system, core::Project(), "", reporter);
        Recording recording;
        std::vector<std::unique_ptr<RecorderDescriptor>> descriptors;
        std::vector<std::unique_ptr<Analysis>> analyses;
        for (const std::string& label : labels) {
            descriptors.push_back(std::make_unique<RecorderDescriptor>(label));
            analyses.push_back(std::make_unique<Recorder>(*descriptors.back(), context, recording));
        }
        const RunResult run = runAnalyses(context, analyses, frames);

        EXPECT_EQ(logs, strings(at(c, "logs")));
        const Array& expectedProgress = list(at(c, "progress"));
        ASSERT_EQ(progress.size(), expectedProgress.size());
        for (std::size_t i = 0; i < progress.size(); ++i) {
            EXPECT_EQ(progress[i].first, hexFloat(list(expectedProgress[i])[0]));
            EXPECT_EQ(progress[i].second, str(list(expectedProgress[i])[1]));
        }
        ASSERT_EQ(run.results.size(), labels.size());
        for (std::size_t i = 0; i < labels.size(); ++i) {
            EXPECT_EQ(run.results[i].name, str(list(at(c, "results"))[i]));
        }
        const Array& times = list(at(c, "times"));
        ASSERT_EQ(run.timePs.size(), times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            EXPECT_EQ(run.timePs[i], hexFloat(times[i])) << i;
        }
        if (const Value* calls = c.find("calls")) {
            std::vector<std::string> expected;
            for (const Value& call : list(*calls)) {
                const Array& parts = list(call);
                std::string text = str(parts[0]) + " " + str(parts[1]) + " ";
                if (str(parts[0]) == "setup") {
                    const auto fs = ints<std::int64_t>(parts[2]);
                    for (std::size_t i = 0; i < fs.size(); ++i) {
                        text += (i ? "," : "") + std::to_string(fs[i]);
                    }
                } else if (str(parts[0]) == "process") {
                    text += std::to_string(integer(parts[2])) + " " + std::to_string(integer(parts[3]));
                } else {
                    text += std::to_string(integer(parts[2]));
                }
                expected.push_back(text);
            }
            EXPECT_EQ(recording.calls, expected);
        } else {
            EXPECT_EQ(static_cast<std::int64_t>(recording.calls.size()), integer(at(c, "n_calls")));
        }
    }
}

TEST(AnalysisFramework, RunAnalysesCancelledAndWithoutTrajectory) {
    core::System system = memorySystem(5, 1.0);
    const core::Reporter reporter;
    reporter.cancel();
    AnalysisContext context(system, core::Project(), "", reporter);
    Recording recording;
    const RecorderDescriptor d("R");
    std::vector<std::unique_ptr<Analysis>> analyses;
    analyses.push_back(std::make_unique<Recorder>(d, context, recording));
    const std::vector<std::int64_t> frames{0, 1};
    EXPECT_THROW((void)runAnalyses(context, analyses, frames), core::Cancelled);
    EXPECT_EQ(recording.calls, (std::vector<std::string>{"setup R 0,1"}));

    core::System bare{makeTopology({0}), nullptr};
    AnalysisContext noFrames(bare, core::Project());
    EXPECT_THROW((void)runAnalyses(noFrames, analyses, frames), std::logic_error);
}

TEST(AnalysisFramework, Context) {
    core::System system = memorySystem(2, 1.0);
    core::Frame frame;
    frame.box = std::array<float, 6>{10, 20, 30, 90, 90, 90};
    core::Project project;
    project.groups = {{"First", "1", "", "#000000"}};
    std::vector<std::string> lines;
    const core::Reporter reporter([&](const std::string& m) { lines.push_back(m); });
    {
        test::LogCapture log;
        const AnalysisContext context(system, project, "top.parm7", reporter);
        EXPECT_EQ(context.box(frame), std::nullopt);  // use_pbc off
        EXPECT_EQ(context.topologyPath(), "top.parm7");
        context.warn("careful");
        // log.warning(message), then the reporter's raw_log (logged at debug level)
        EXPECT_EQ(log.lines(), (std::vector<std::string>{"WARNING: careful", "DEBUG: WARNING: careful"}));
        EXPECT_EQ(lines, (std::vector<std::string>{"WARNING: careful"}));
        EXPECT_EQ(context.groupAtoms("First").indices(), (std::vector<std::int32_t>{0}));
        EXPECT_EQ(context.select("name A0").size(), 1u);
        EXPECT_EQ(context.select("name Z", true).size(), 0u);
        EXPECT_THROW((void)context.select("name Z"), core::ValueError);
        EXPECT_THROW((void)context.groupAtoms("Nope"), core::ValueError);
    }
    project.analyses.use_pbc = true;
    const AnalysisContext periodic(system, project);
    EXPECT_EQ(periodic.box(frame), frame.box);
    EXPECT_TRUE(periodic.usePbc());
    EXPECT_EQ(periodic.box(core::Frame()), std::nullopt);
}
