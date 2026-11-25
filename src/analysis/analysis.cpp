#include "mdbridge/analysis/analysis.hpp"

#include <stdexcept>

#include "builtin.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/hash.hpp"

namespace mdbridge::analysis {

namespace detail {

core::Table sectionParams(const core::Project& project, std::string_view section) {
    const core::Value plain = core::toPlain(project);
    const core::Value* entry = plain.path({"analyses", section});
    const core::Table* table = entry ? entry->get<core::Table>() : nullptr;
    if (!table) {
        throw std::logic_error("no analyses." + std::string(section) + " section");
    }
    core::Table out = *table;
    for (auto& [key, value] : out) {
        if (key == "enabled") {
            value = core::Value();
        }
    }
    return out;
}

GroupFields groupFields(const core::Project& project, std::string_view section) {
    const core::AnalysisParams& a = project.analyses;
    if (section == "contacts") {
        return {&a.contacts.group_a, &a.contacts.group_b};
    }
    if (section == "hbonds") {
        return {&a.hbonds.group_a, &a.hbonds.group_b};
    }
    if (section == "mg") {
        return {&a.mg.group_a, &a.mg.group_b};
    }
    if (section == "stacking") {
        return {&a.stacking.group_a, &a.stacking.group_b};
    }
    if (section == "energy") {
        return {&a.energy.group_a, &a.energy.group_b};
    }
    return {};
}

MutableGroupFields groupFields(core::Project& project, std::string_view section) {
    const GroupFields fields = groupFields(static_cast<const core::Project&>(project), section);
    return {const_cast<std::string*>(fields.a), const_cast<std::string*>(fields.b)};
}

void notImplemented(const AnalysisDescriptor& descriptor) {
    throw core::NamedError("NotImplementedError", descriptor.name() + ": not implemented yet");
}

}  // namespace detail

// ---------------------------------------------------------------- ExpectedOutputs

const std::vector<std::string>& ExpectedOutputs::of(std::string_view kind) const {
    if (kind == "series") {
        return series;
    }
    if (kind == "matrices") {
        return matrices;
    }
    if (kind == "pairmaps") {
        return pairmaps;
    }
    throw std::out_of_range(std::string(kind));
}

// ---------------------------------------------------------------- AnalysisDescriptor

AnalysisDescriptor::AnalysisDescriptor(Info info) : m_info(std::move(info)) {}

bool AnalysisDescriptor::enabled(const core::Project&) const {
    return false;
}

core::Table AnalysisDescriptor::paramsFor(const core::Project&) const {
    return {};
}

std::vector<std::string> AnalysisDescriptor::selectionTexts(const core::Project&) const {
    return {};
}

ExpectedOutputs AnalysisDescriptor::expectedOutputs(const core::Project&) const {
    return {};
}

std::vector<std::string> AnalysisDescriptor::groups(const core::Project& project) const {
    const detail::GroupFields fields = detail::groupFields(project, name());
    if (groupSlots() == 2 && fields.a) {
        return {*fields.a, *fields.b};
    }
    return {};
}

void AnalysisDescriptor::withGroups(core::Project& project, const std::vector<std::string>& groups) const {
    const detail::MutableGroupFields fields = detail::groupFields(project, name());
    if (!fields.a) {
        throw std::logic_error("analyses." + name() + " has no residue groups");
    }
    if (groups.size() != 2) {
        throw std::logic_error(name() + " takes 2 groups");
    }
    *fields.a = groups[0];
    *fields.b = groups[1];
}

std::string AnalysisDescriptor::outputLabel(const core::Project& project, std::string_view out) const {
    std::vector<std::string> names = groups(project);
    names.emplace_back("A");
    names.emplace_back("B");
    for (const auto& [key, text] : outputLabels()) {
        if (key != out) {
            continue;
        }
        // str.format(a=..., b=...) of a template whose only fields are {a} and {b}
        std::string filled;
        for (std::size_t pos = 0; pos < text.size(); ++pos) {
            if (text.compare(pos, 3, "{a}") == 0) {
                filled += names[0];
                pos += 2;
            } else if (text.compare(pos, 3, "{b}") == 0) {
                filled += names[1];
                pos += 2;
            } else {
                filled += text[pos];
            }
        }
        return filled;
    }
    return label() + ": " + std::string(out);
}

std::string AnalysisDescriptor::cacheKey(const core::Project& project, const core::Value& inputs) const {
    const core::Table key{{"analysis", name()},
                          {"version", version()},
                          {"params", paramsFor(project)},
                          {"groups", referencedGroups(selectionTexts(project), project)},
                          {"inputs", inputs},
                          {"use_pbc", project.analyses.use_pbc}};
    return core::stableHash(key);
}

// ---------------------------------------------------------------- Analysis

Analysis::Analysis(const AnalysisDescriptor& descriptor, AnalysisContext& context, std::string instance)
    : m_descriptor(&descriptor), m_context(&context),
      m_instance(instance.empty() ? descriptor.name() : std::move(instance)) {}

void Analysis::prepare(const core::Frame&) {}

void Analysis::setup(std::span<const std::int64_t> frames) {
    m_frameCount = frames.size();
}

core::AnalysisResult Analysis::newResult(const std::string& key, core::Table meta) const {
    const auto assign = [&meta](const std::string& k, core::Value v) {
        for (auto& [entry, value] : meta) {
            if (entry == k) {
                value = std::move(v);
                return;
            }
        }
        meta.emplace_back(k, std::move(v));
    };
    assign("label", label());
    assign("analysis", name());
    core::AnalysisResult result;
    result.name = m_instance;
    result.key = key;
    result.params = m_descriptor->paramsFor(project());
    result.meta = std::move(meta);
    return result;
}

}  // namespace mdbridge::analysis
