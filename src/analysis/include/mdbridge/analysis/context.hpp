#pragma once
// Shared state handed to every analysis of one run (mdatk/analysis/base.py AnalysisContext): the
// loaded system, the project the analysis runs with (an instance's own copy, see
// registry.hpp instanceProject), the topology file and the reporter.

#include <array>
#include <optional>
#include <string>
#include <string_view>

#include "mdbridge/core/progress.hpp"
#include "mdbridge/core/project.hpp"
#include "mdbridge/core/topology.hpp"

namespace mdbridge::analysis {

class AnalysisContext {
public:
    // `system` is shared by the contexts of one run and must outlive this context.
    AnalysisContext(core::System& system, core::Project project, std::string topologyPath = {},
                    core::Reporter reporter = core::Reporter());

    core::System& system() const { return *m_system; }
    const core::Topology& topology() const { return *m_system->topology; }
    const core::Project& project() const { return m_project; }
    const std::string& topologyPath() const { return m_topologyPath; }
    const core::Reporter& reporter() const { return m_reporter; }
    // analyses.use_pbc of the project
    bool usePbc() const { return m_usePbc; }

    // select(u, text, project.groups, allow_empty) on the system's current frame (runAnalyses makes
    // that the first analysed frame for every setup, SEL-17): sorted unique atoms; distance
    // keywords use the frame's cell only with analyses.use_pbc (SEL-17; MDATK: always). Throws
    // core::SelectionError.
    core::AtomGroup select(std::string_view text, bool allowEmpty = false) const;
    // The atoms of group:<name>, heavy atoms only when asked.
    core::AtomGroup groupAtoms(std::string_view name, bool heavyOnly = false) const;
    // The frame's unit cell when analyses.use_pbc is set (std::nullopt otherwise or without a cell).
    std::optional<std::array<float, 6>> box(const core::Frame& frame) const;
    // A warning: logged (core::logWarning) and passed to the reporter as "WARNING: <message>".
    void warn(const std::string& message) const;

private:
    core::System* m_system;
    core::Project m_project;
    std::string m_topologyPath;
    core::Reporter m_reporter;
    bool m_usePbc;
};

}  // namespace mdbridge::analysis
