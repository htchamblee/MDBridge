#include "mdbridge/analysis/context.hpp"

#include <utility>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/selections.hpp"

namespace mdbridge::analysis {

AnalysisContext::AnalysisContext(core::System& system, core::Project project, std::string topologyPath,
                                 core::Reporter reporter)
    : m_system(&system), m_project(std::move(project)), m_topologyPath(std::move(topologyPath)),
      m_reporter(std::move(reporter)), m_usePbc(m_project.analyses.use_pbc) {}

core::AtomGroup AnalysisContext::select(std::string_view text, bool allowEmpty) const {
    return core::select(*m_system, text, m_project.groups, allowEmpty);
}

core::AtomGroup AnalysisContext::groupAtoms(std::string_view name, bool heavyOnly) const {
    core::AtomGroup atoms = select("group:" + std::string(name));
    return heavyOnly ? core::heavy(atoms) : atoms;
}

std::optional<std::array<float, 6>> AnalysisContext::box(const core::Frame& frame) const {
    if (!m_usePbc) {
        return std::nullopt;
    }
    return frame.box;
}

void AnalysisContext::warn(const std::string& message) const {
    core::logWarning(message);
    m_reporter.rawLog("WARNING: " + message);
}

}  // namespace mdbridge::analysis
