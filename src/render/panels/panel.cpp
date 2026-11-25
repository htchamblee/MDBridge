// The panel base (panel.hpp): mdatk/render/panels.py RenderContext and Panel (MDATK's own code, ported).

#include "mdbridge/render/panels/panel.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"

namespace mdbridge::render::panels {

Context::Context(const core::Dataset& ds, const core::Project& p)
    : dataset(ds),
      project(p),
      theme(render::theme(p.visual.theme)),
      time(ds.time(p.inputs.time_unit)),
      timeUnit(p.inputs.time_unit),
      timeLabel(ds.timeLabelFor(p.inputs.time_unit)),
      fontSize(p.visual.font_size) {}

std::string Context::formatTime(std::size_t pos) const {
    const std::int64_t frame = dataset.frames().at(pos);
    if (timeUnit == "frame") return "frame " + std::to_string(frame);
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.3f", time.at(pos));
    std::string unit = timeUnit;
    if (const std::size_t at = unit.find("us"); at != std::string::npos) unit.replace(at, 2, "\xc2\xb5s");
    return "t = " + std::string(buf) + " " + unit + "  (frame " + std::to_string(frame) + ")";
}

namespace {

// Dataset lookups throw std::out_of_range with Python's KeyError text
template <class F>
auto lookup(F&& f) {
    try {
        return f();
    } catch (const std::out_of_range& e) {
        throw NoData(e.what());
    }
}

}  // namespace

core::Series Context::series(std::string_view key) const {
    return lookup([&] { return dataset.series(key); });
}
core::Matrix Context::matrix(std::string_view key) const {
    return lookup([&] { return dataset.matrix(key); });
}
core::PairMap Context::pairmap(std::string_view key) const {
    return lookup([&] { return dataset.pairmap(key); });
}

std::string ArtistRef::id() const {
    return "ax" + std::to_string(axes) + "/" + (legend ? std::string("legend") : std::to_string(index));
}

Panel::Panel(std::string type, const core::PanelSpec& spec, Context& context)
    : ctx(context), theme(context.theme), m_type(std::move(type)), m_spec(spec) {}

void Panel::style(Axes& a, bool grid) const { styleAxes(a, axesStyle(theme, grid, ctx.fontSize)); }

void Panel::init(Figure& fig, SubplotSpec slot) {
    m_axes = fig.addSubplot(slot);
    style(fig.axes[m_axes]);
    try {
        setup(fig);
    } catch (const NoData& e) {
        showMessage(fig, std::string("No data: ") + e.what());
    }
    const std::string title = !m_spec.title.empty() ? m_spec.title : m_autoTitle;
    if (!title.empty()) {
        Text& t = fig.axes[m_axes].setTitle(title, HAlign::Left);
        t.color = toRgba(theme.ink);
        t.font.size = ctx.fontSize;
    }
}

void Panel::showMessage(Figure& fig, const std::string& text) {
    m_ok = false;
    // what setup added beyond the main axes (a twin or colour bar made before the KeyError) stays,
    // as in MDATK; the main axes is cleared
    Axes& a = fig.axes[m_axes];
    a.clear();
    a.face = toRgba(theme.surface);
    a.setXTicks({});
    a.setYTicks({});
    for (Spine& s : a.spines) s.color = toRgba(theme.grid);
    Text t;
    t.text = text;
    t.xy = {0.5, 0.5};
    t.coords = Coords::Axes;
    t.ha = HAlign::Center;
    t.va = VAlign::Center;
    t.wrap = true;
    t.color = toRgba(theme.muted);
    t.font.size = ctx.fontSize * 0.85;
    a.add(std::move(t));
}

const core::Value* Panel::option(std::string_view key) const {
    for (const auto& [k, v] : m_spec.options) {
        if (k == key) return &v;
    }
    return nullptr;
}

namespace {

[[noreturn]] void badOption(std::string_view key, std::string_view expected, const core::Value& v) {
    throw core::ValueError("panel option '" + std::string(key) + "': expected " + std::string(expected) + ", got " +
                           core::pythonRepr(v));
}

}  // namespace

double Panel::optionNumber(std::string_view key, double fallback) const {
    const core::Value* v = option(key);
    if (!v || v->isNull()) return fallback;
    if (const auto* d = v->get<double>()) return *d;
    if (const auto* i = v->get<std::int64_t>()) return static_cast<double>(*i);
    if (const auto* b = v->get<bool>()) return *b ? 1.0 : 0.0;
    if (const auto* s = v->get<std::string>()) {
        if (const auto x = core::tryPythonFloat(*s)) return *x;
    }
    badOption(key, "a number", *v);
}

std::int64_t Panel::optionInt(std::string_view key, std::int64_t fallback) const {
    const core::Value* v = option(key);
    if (!v || v->isNull()) return fallback;
    if (const auto* i = v->get<std::int64_t>()) return *i;
    if (const auto* b = v->get<bool>()) return *b ? 1 : 0;
    if (const auto* d = v->get<double>()) {
        if (std::isfinite(*d)) return static_cast<std::int64_t>(std::trunc(*d));  // int(5.7) == 5
    }
    if (const auto* s = v->get<std::string>()) {
        if (const auto x = core::tryPythonInt(*s)) return *x;
    }
    badOption(key, "an integer", *v);
}

bool Panel::optionBool(std::string_view key, bool fallback) const {
    const core::Value* v = option(key);
    if (!v) return fallback;
    // Python's bool(): truthiness
    if (v->isNull()) return false;
    if (const auto* b = v->get<bool>()) return *b;
    if (const auto* i = v->get<std::int64_t>()) return *i != 0;
    if (const auto* d = v->get<double>()) return *d != 0.0;
    if (const auto* s = v->get<std::string>()) return !s->empty();
    if (v->isArray() || v->isTable()) return v->size() != 0;
    return true;
}

std::string Panel::optionString(std::string_view key, std::string fallback) const {
    const core::Value* v = option(key);
    if (!v) return fallback;
    if (const auto* s = v->get<std::string>()) return *s;
    return core::pythonStr(*v);
}

}  // namespace mdbridge::render::panels
