// The rendering part of check() (checks.hpp; RN-11 ... RN-15).

#include "mdbridge/render/panels/checks.hpp"

#include <cmath>
#include <cstdint>
#include <string_view>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/render/color.hpp"
#include "mdbridge/render/colormap.hpp"
#include "mdbridge/render/panels/charts.hpp"

namespace mdbridge::render::panels {

namespace {

const core::Value* option(const core::Table& options, std::string_view key) {
    for (const auto& [k, v] : options) {
        if (k == key) return &v;
    }
    return nullptr;
}

// what Panel::optionNumber accepts: null, numbers, booleans, text Python's float() reads
bool readsAsNumber(const core::Value& v) {
    if (v.isNull() || v.get<double>() || v.get<std::int64_t>() || v.get<bool>()) return true;
    if (const auto* s = v.get<std::string>()) return core::tryPythonFloat(*s).has_value();
    return false;
}

// what Panel::optionInt accepts: null, integers, booleans, finite floats (truncated), text int() reads
std::optional<std::int64_t> readsAsInt(const core::Value& v, bool& ok) {
    ok = true;
    if (v.isNull()) return std::nullopt;
    if (const auto* i = v.get<std::int64_t>()) return *i;
    if (const auto* b = v.get<bool>()) return *b ? 1 : 0;
    if (const auto* d = v.get<double>()) {
        if (std::isfinite(*d)) return static_cast<std::int64_t>(std::trunc(*d));
    }
    if (const auto* s = v.get<std::string>()) {
        if (const auto x = core::tryPythonInt(*s)) return *x;
    }
    ok = false;
    return std::nullopt;
}

class PanelChecker {
public:
    PanelChecker(std::vector<std::string>& out, std::size_t k, const core::Table& options)
        : m_out(out), m_head("layout.panels[" + std::to_string(k) + "].options."), m_options(options) {}

    void number(std::string_view key) {
        if (const core::Value* v = option(m_options, key); v && !readsAsNumber(*v)) {
            problem(key, "expected a number, got " + core::pythonRepr(*v));
        }
    }

    void integer(std::string_view key, std::optional<std::int64_t> minimum = std::nullopt) {
        const core::Value* v = option(m_options, key);
        if (!v) return;
        bool ok = false;
        const std::optional<std::int64_t> value = readsAsInt(*v, ok);
        if (!ok) {
            problem(key, "expected a whole number, got " + core::pythonRepr(*v));
        } else if (minimum && value && *value < *minimum) {
            problem(key, "expected a whole number of at least " + std::to_string(*minimum) + ", got " +
                             core::pythonRepr(*v));
        }
    }

    void pair(std::string_view key) {
        const core::Value* v = option(m_options, key);
        if (!v) return;
        const core::Array* a = v->get<core::Array>();
        bool ok = a && a->size() == 2;
        if (ok) {
            for (const core::Value& x : *a) {
                ok = ok && (x.get<double>() || x.get<std::int64_t>());
            }
        }
        if (!ok) problem(key, "expected two numbers, got " + core::pythonRepr(*v));
    }

    void choice(std::string_view key, std::initializer_list<std::string_view> choices) {
        const core::Value* v = option(m_options, key);
        if (!v) return;
        const std::string* s = v->get<std::string>();
        for (std::string_view c : choices) {
            if (s && *s == c) return;
        }
        std::string list;
        std::size_t i = 0;
        for (std::string_view c : choices) {
            list += (i == 0 ? "" : i + 1 == choices.size() ? " or " : ", ") + std::string(c);
            ++i;
        }
        problem(key, "expected " + list + ", got " + core::pythonRepr(*v));
    }

private:
    void problem(std::string_view key, const std::string& text) { m_out.push_back(m_head + std::string(key) + ": " + text); }

    std::vector<std::string>& m_out;
    std::string m_head;
    const core::Table& m_options;
};

void color(std::vector<std::string>& out, const std::string& field, const std::string& value) {
    if (const auto problem = colorProblem(value)) out.push_back(field + ": " + *problem);
}

}  // namespace

std::vector<std::string> renderProblems(const core::Project& project) {
    std::vector<std::string> out;
    // RN-11
    if (const auto problem = charts::chartFormatProblem(project.output.charts.formats)) {
        out.push_back("output.charts.formats: " + *problem);
    }
    // RN-12
    for (std::size_t k = 0; k < project.groups.size(); ++k) {
        color(out, "groups[" + std::to_string(k) + "].color", project.groups[k].color);
    }
    const core::VisualParams& v = project.visual;
    color(out, "visual.base_color", v.base_color);
    color(out, "visual.proximity_color", v.proximity_color);
    color(out, "visual.glow_color", v.glow_color);
    color(out, "visual.ion_color", v.ion_color);
    if (!project.density.map_path.empty()) color(out, "density.color", project.density.color);
    if (!project.output.charts.color.empty()) color(out, "output.charts.color", project.output.charts.color);
    // RN-13
    if (const auto problem = colormapProblem(v.chain_palette)) out.push_back("visual.chain_palette: " + *problem);
    // RN-15 (RN-04: the heatmap's limits)
    const std::vector<core::PanelSpec>& panels = project.layout.panels;
    for (std::size_t k = 0; k < panels.size(); ++k) {
        PanelChecker check(out, k, panels[k].options);
        const std::string& type = panels[k].type;
        if (type == "timeseries") {
            check.integer("smooth");
            check.number("target");
            check.pair("ylim");
        } else if (type == "heatmap") {
            check.number("vmin");
            check.number("vmax");
            check.integer("max_ticks", 1);
        } else if (type == "contact_map") {
            check.choice("mode", {"instant", "frame", "frequency"});
            check.integer("frame");
            check.integer("max_ticks", 1);
        }
    }
    return out;
}

std::string figureFileFormat(const std::string& path) {
    const std::string ext = core::posixpath::splitext(core::posixpath::basename(path)).second;
    return ext.empty() ? std::string("png") : charts::normalizeFormat(ext);
}

std::optional<std::string> figureFileProblem(const std::string& path, const std::string& what) {
    const std::string format = figureFileFormat(path);
    std::string list;
    for (std::string_view f : charts::chartFormats()) {
        if (f == format) return std::nullopt;
        list += (list.empty() ? "" : ", ") + std::string(f);
    }
    return what + " format " + core::pythonStringRepr(format) + " is not supported (choose from " + list + ")";
}

}  // namespace mdbridge::render::panels
