#pragma once
// The analysis framework (mdatk/analysis/base.py): the static description of an analysis type
// (name, label, parameters that feed the cache key, outputs, residue-group slots), the per-run
// streaming analysis object (setup / process every frame / finalize) and the single-pass runner.
//
// To add an analysis: derive a descriptor from AnalysisDescriptor and an analysis from Analysis,
// and list the descriptor in registry.cpp (the registry order is MDATK's import order).

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/analysis/context.hpp"
#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"
#include "mdbridge/core/value.hpp"

namespace mdbridge::analysis {

class Analysis;

// Output names of an analysis per kind, before running ({"series", "matrices", "pairmaps"}).
struct ExpectedOutputs {
    std::vector<std::string> series;
    std::vector<std::string> matrices;
    std::vector<std::string> pairmaps;

    // The list of a kind name ("series", "matrices", "pairmaps"); std::out_of_range otherwise.
    const std::vector<std::string>& of(std::string_view kind) const;
    bool operator==(const ExpectedOutputs&) const = default;
};

// Output name -> label template; "{a}" / "{b}" are replaced by the group names.
using OutputLabels = std::vector<std::pair<std::string, std::string>>;

// The class-level part of an MDATK analysis (its attributes and classmethods). One static instance
// per analysis type, owned by the registry.
class AnalysisDescriptor {
public:
    struct Info {
        std::string name;
        std::string label;
        int version = 1;
        bool internal = false;  // trace / proximity / atoms / basepairs feed panels only
        int groupSlots = 0;  // residue groups an instance may choose (0, 1 or 2)
        OutputLabels outputLabels;
    };

    explicit AnalysisDescriptor(Info info);
    virtual ~AnalysisDescriptor() = default;
    AnalysisDescriptor(const AnalysisDescriptor&) = delete;
    AnalysisDescriptor& operator=(const AnalysisDescriptor&) = delete;

    const std::string& name() const { return m_info.name; }
    const std::string& label() const { return m_info.label; }
    int version() const { return m_info.version; }
    bool internal() const { return m_info.internal; }
    int groupSlots() const { return m_info.groupSlots; }
    const OutputLabels& outputLabels() const { return m_info.outputLabels; }

    // enabled(project); false by default (internal analyses).
    virtual bool enabled(const core::Project& project) const;
    // params_for(project): the parameters that influence the result (cache key, result params).
    // Entries MDATK sets to None (e.g. "enabled") are Null here; the cache key drops them.
    virtual core::Table paramsFor(const core::Project& project) const;
    // selection_texts(project): selections and group names used (for the cache key's groups).
    virtual std::vector<std::string> selectionTexts(const core::Project& project) const;
    // expected_outputs(project)
    virtual ExpectedOutputs expectedOutputs(const core::Project& project) const;
    // groups(project): (group_a, group_b) of analyses.<name> for two slots, else none.
    virtual std::vector<std::string> groups(const core::Project& project) const;
    // with_groups(project, groups): sets analyses.<name>.group_a / group_b (two slots); `project`
    // is the caller's private copy.
    virtual void withGroups(core::Project& project, const std::vector<std::string>& groups) const;
    // output_label(project, out): the output's template with {a} / {b} filled from groups() (or
    // "A" / "B"), else "<label>: <out>".
    virtual std::string outputLabel(const core::Project& project, std::string_view out) const;

    // cache_key(project, inputs): stableHash({"analysis", "version", "params",
    // "groups": referencedGroups(selectionTexts), "inputs", "use_pbc"}), 16 hex digits.
    std::string cacheKey(const core::Project& project, const core::Value& inputs) const;

    // The streaming analysis for one run (`instance`: the instance id, "" for the name).
    virtual std::unique_ptr<Analysis> create(AnalysisContext& context, std::string instance = {}) const = 0;

private:
    Info m_info;
};

// The per-run part of an MDATK analysis: frames are streamed once, in order.
class Analysis {
public:
    Analysis(const AnalysisDescriptor& descriptor, AnalysisContext& context, std::string instance = {});
    virtual ~Analysis() = default;
    Analysis(const Analysis&) = delete;
    Analysis& operator=(const Analysis&) = delete;

    const AnalysisDescriptor& descriptor() const { return *m_descriptor; }
    const std::string& name() const { return m_descriptor->name(); }
    const std::string& label() const { return m_descriptor->label(); }
    // The instance id ("contacts@A+B"), the name when none was given.
    const std::string& instance() const { return m_instance; }
    AnalysisContext& context() const { return *m_context; }
    const core::Project& project() const { return m_context->project(); }
    // Frames given to setup().
    std::size_t frameCount() const { return m_frameCount; }

    // Called once with the analysed frame numbers before the first process(); the base stores
    // their count.
    virtual void setup(std::span<const std::int64_t> frames);
    // Frame k (0-based position among the analysed frames).
    virtual void process(std::size_t k, const core::Frame& frame) = 0;
    virtual core::AnalysisResult finalize() = 0;

    // new_result(key, **meta): name = instance, params = paramsFor(project), meta = meta with
    // "label" and "analysis" set (replaced in place when present, else appended).
    core::AnalysisResult newResult(const std::string& key = {}, core::Table meta = {}) const;

private:
    const AnalysisDescriptor* m_descriptor;
    AnalysisContext* m_context;
    std::string m_instance;
    std::size_t m_frameCount = 0;
};

struct RunResult {
    std::vector<core::AnalysisResult> results;  // in the order of the analyses
    std::vector<double> timePs;  // time of every analysed frame
};

// run_analyses(ctx, analyses, frame_indices, reporter): logs "Preparing <label> ..." and calls
// setup() for each analysis, logs "Analysing N frame(s): <labels joined by ', '>", then seeks
// every frame of `frames` in order on the context's trajectory: checkCancel(), the frame's time,
// process(k, frame) of each analysis in order, and progress((k + 1) / n, "Analysing frame k+1/n")
// when k % max(1, n / 200) == 0 or on the last frame. Then finalize() of each. `reporter`:
// the context's when null. The trajectory is left on the last analysed frame. Throws
// std::logic_error when the system has no trajectory, core::Cancelled when cancelled, and
// whatever the analyses or the reader throw.
RunResult runAnalyses(AnalysisContext& context, std::span<Analysis* const> analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter = nullptr);
RunResult runAnalyses(AnalysisContext& context, const std::vector<std::unique_ptr<Analysis>>& analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter = nullptr);

}  // namespace mdbridge::analysis
