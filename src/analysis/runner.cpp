#include <algorithm>
#include <stdexcept>
#include <string>

#include "mdbridge/analysis/analysis.hpp"

namespace mdbridge::analysis {

RunResult runAnalyses(AnalysisContext& context, std::span<Analysis* const> analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter) {
    const core::Reporter& rep = reporter ? *reporter : context.reporter();
    core::Trajectory* trajectory = context.system().trajectory.get();
    if (!trajectory) {
        throw std::logic_error("runAnalyses: the system has no trajectory");
    }
    const std::size_t n = frames.size();
    std::string labels;  // ', '.join(a.label for a in analyses)
    for (std::size_t i = 0; i < analyses.size(); ++i) {
        rep.log("Preparing " + analyses[i]->label() + " ...");
        analyses[i]->setup(frames);
        labels += (i ? ", " : "") + analyses[i]->label();
    }
    RunResult out;
    out.timePs.resize(n);
    const std::size_t every = std::max<std::size_t>(1, n / 200);
    rep.log("Analysing " + std::to_string(n) + " frame(s): " + labels);
    for (std::size_t k = 0; k < n; ++k) {
        rep.checkCancel();
        const core::Frame& frame = trajectory->seek(frames[k]);
        out.timePs[k] = frame.time;
        for (Analysis* a : analyses) {
            a->process(k, frame);
        }
        if (k % every == 0 || k == n - 1) {
            rep.progress(static_cast<double>(k + 1) / static_cast<double>(n),
                         "Analysing frame " + std::to_string(k + 1) + "/" + std::to_string(n));
        }
    }
    out.results.reserve(analyses.size());
    for (Analysis* a : analyses) {
        out.results.push_back(a->finalize());
    }
    return out;
}

RunResult runAnalyses(AnalysisContext& context, const std::vector<std::unique_ptr<Analysis>>& analyses,
                      std::span<const std::int64_t> frames, const core::Reporter* reporter) {
    std::vector<Analysis*> pointers;
    pointers.reserve(analyses.size());
    for (const auto& a : analyses) {
        pointers.push_back(a.get());
    }
    return runAnalyses(context, pointers, frames, reporter);
}

}  // namespace mdbridge::analysis
