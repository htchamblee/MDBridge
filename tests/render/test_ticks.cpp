// Axis limits, tick locations and tick labels (src/render/core/ticks.cpp) against what Matplotlib
// 3.11.2 and MDATK compute (data/render/ticks/tick_cases.json, tests/parity/make_tick_cases.py):
// every number bit for bit (the sign of zero included), every string exactly.

#include <cmath>
#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/render/ticks.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using core::Array;
using core::Value;
using render::AxisDirection;
using render::MaxNLocator;
using render::MaxNLocatorParams;
using render::Range;
using render::ScalarFormatter;
using render::ScalarFormatterOptions;

namespace {

const Value& data() {
    static const Value value = test::readJson(test::kData / "render" / "ticks" / "tick_cases.json");
    return value;
}

const Value& at(const Value& v, std::string_view key) {
    const Value* found = v.find(key);
    if (!found) {
        throw std::runtime_error("missing key " + std::string(key));
    }
    return *found;
}
bool has(const Value& v, std::string_view key) { return v.find(key) != nullptr && !v.find(key)->isNull(); }

double num(const Value& v) {
    if (const auto* d = v.get<double>()) {
        return *d;
    }
    if (const auto* i = v.get<std::int64_t>()) {
        return static_cast<double>(*i);
    }
    throw std::runtime_error("not a number");
}
long long integer(const Value& v) { return *v.get<std::int64_t>(); }
bool flag(const Value& v) { return *v.get<bool>(); }
std::string str(const Value& v) { return *v.get<std::string>(); }
const Array& list(const Value& v) { return *v.get<Array>(); }

std::vector<double> nums(const Value& v) {
    std::vector<double> out;
    for (const Value& x : list(v)) {
        out.push_back(num(x));
    }
    return out;
}
std::vector<std::string> strs(const Value& v) {
    std::vector<std::string> out;
    for (const Value& x : list(v)) {
        out.push_back(str(x));
    }
    return out;
}
Range range(const Value& v) {
    const auto xs = nums(v);
    return {xs.at(0), xs.at(1)};
}

bool same(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::isnan(a) && std::isnan(b);
    }
    return a == b && std::signbit(a) == std::signbit(b);
}
bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!same(a[i], b[i])) {
            return false;
        }
    }
    return true;
}
bool same(Range a, Range b) { return same(a.lo, b.lo) && same(a.hi, b.hi); }

std::string show(const std::vector<double>& xs) {
    std::ostringstream out;
    out.precision(17);
    out << '[';
    for (std::size_t i = 0; i < xs.size(); ++i) {
        out << (i ? ", " : "") << xs[i];
    }
    out << ']';
    return out.str();
}
std::string show(Range r) { return show(std::vector<double>{r.lo, r.hi}); }
std::string show(const std::vector<std::string>& xs) {
    std::string out = "[";
    for (std::size_t i = 0; i < xs.size(); ++i) {
        out += (i ? ", '" : "'") + xs[i] + "'";
    }
    return out + "]";
}

std::string caseText(const Value& c) { return core::toJson(c, -1).substr(0, 400); }

MaxNLocator locatorFrom(const Value& p) {
    if (const Value* kind = p.find("kind")) {
        return str(*kind) == "auto" ? MaxNLocator::autoLocator() : MaxNLocator::integerLocator();
    }
    MaxNLocatorParams params;
    const Value& nbins = at(p, "nbins");
    if (nbins.isString()) {
        params.nbins = std::nullopt;
    } else {
        params.nbins = static_cast<int>(integer(nbins));
    }
    if (has(p, "steps")) {
        params.steps = nums(at(p, "steps"));
    }
    if (has(p, "integer")) {
        params.integer = flag(at(p, "integer"));
    }
    if (has(p, "symmetric")) {
        params.symmetric = flag(at(p, "symmetric"));
    }
    if (has(p, "min_n_ticks")) {
        params.minNTicks = static_cast<int>(integer(at(p, "min_n_ticks")));
    }
    if (has(p, "prune")) {
        const std::string prune = str(at(p, "prune"));
        params.prune = prune == "lower" ? render::Prune::Lower : prune == "upper" ? render::Prune::Upper
                                                                                 : render::Prune::Both;
    }
    return MaxNLocator(std::move(params));
}

ScalarFormatterOptions formatterOptions(const Value& o) {
    ScalarFormatterOptions options;
    if (const Value* off = o.find("useOffset")) {
        if (const auto* b = off->get<bool>()) {
            options.useOffset = *b;
        } else {
            options.fixedOffset = num(*off);
        }
    }
    if (has(o, "useMathText")) {
        options.useMathText = flag(at(o, "useMathText"));
    }
    if (has(o, "unicode_minus")) {
        options.unicodeMinus = flag(at(o, "unicode_minus"));
    }
    if (has(o, "powerlimits")) {
        const auto lims = nums(at(o, "powerlimits"));
        options.powerLimitLow = static_cast<int>(lims[0]);
        options.powerLimitHigh = static_cast<int>(lims[1]);
    }
    if (has(o, "scientific")) {
        options.scientific = flag(at(o, "scientific"));
    }
    if (has(o, "offset_threshold")) {
        options.offsetThreshold = static_cast<int>(integer(at(o, "offset_threshold")));
    }
    return options;
}

// Reports the first failures of a section and counts all of them.
class Failures {
public:
    explicit Failures(std::string section) : m_section(std::move(section)) {}
    ~Failures() { EXPECT_EQ(m_count, 0) << m_count << " failing " << m_section << " cases"; }
    void add(std::size_t index, const Value& c, const std::string& what) {
        if (m_count++ < 12) {
            ADD_FAILURE() << m_section << " case " << index << ": " << what << "\n  " << caseText(c);
        }
    }

private:
    std::string m_section;
    int m_count = 0;
};

// One axis of a Figure's Axes as Axis._update_ticks sees it, given its view interval.
void checkAxis(Failures& failures, std::size_t index, const Value& c, const Value& rec, AxisDirection dir,
               double labelSize, Range view) {
    const int space = render::tickSpace(dir, num(at(rec, "length_points")), labelSize);
    if (space != integer(at(rec, "tick_space"))) {
        failures.add(index, c, "tick space " + std::to_string(space));
    }
    if (!same(view, range(at(rec, "view")))) {
        failures.add(index, c, "view " + show(view));
    }
    const MaxNLocator locator = flag(at(rec, "integer")) ? MaxNLocator::integerLocator() : MaxNLocator::autoLocator();
    const std::vector<double> locs = locator.tickValues(view.lo, view.hi, space);
    if (!same(locs, nums(at(rec, "locs")))) {
        failures.add(index, c, "locs " + show(locs));
        return;
    }
    ScalarFormatter formatter;
    const std::vector<std::string> labels = formatter.formatTicks(locs, view);
    if (labels != strs(at(rec, "labels")) || formatter.offsetText() != str(at(rec, "offset_text"))) {
        failures.add(index, c, "labels " + show(labels) + " offset '" + formatter.offsetText() + "'");
    }
    std::vector<double> drawn;
    for (const std::size_t i : render::ticksToDraw(locs, view)) {
        drawn.push_back(static_cast<double>(i));
    }
    if (!same(drawn, nums(at(rec, "drawn")))) {
        failures.add(index, c, "drawn " + show(drawn));
    }
}

}  // namespace

TEST(Ticks, Nonsingular) {
    Failures failures("nonsingular");
    const Array& cases = list(at(data(), "nonsingular"));
    ASSERT_GT(cases.size(), 100u);
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const Range r = render::nonsingular(num(at(c, "vmin")), num(at(c, "vmax")), num(at(c, "expander")),
                                            num(at(c, "tiny")), flag(at(c, "increasing")));
        if (!same(r, range(at(c, "result")))) {
            failures.add(i, c, show(r));
        }
    }
}

TEST(Ticks, MaxNLocator) {
    Failures failures("locator");
    const Array& cases = list(at(data(), "locator"));
    ASSERT_GT(cases.size(), 1000u);
    std::size_t checked = 0;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const Value& p = at(c, "params");
        if (const Value* invalid = c.find("invalid_steps")) {
            if (nums(at(p, "steps")).empty()) {
                continue;  // Matplotlib: IndexError; MDBridge: empty means the default steps
            }
            try {
                (void)locatorFrom(p);
                failures.add(i, c, "accepted");
            } catch (const std::invalid_argument& e) {
                if (e.what() != str(*invalid)) {
                    failures.add(i, c, e.what());
                }
            }
            continue;
        }
        if (has(c, "error")) {
            failures.add(i, c, "unexpected Matplotlib error in the expectations");
            continue;
        }
        const MaxNLocator locator = locatorFrom(p);
        std::optional<int> space;
        if (has(c, "tick_space")) {
            space = static_cast<int>(integer(at(c, "tick_space")));
        }
        const double vmin = num(at(c, "vmin"));
        const double vmax = num(at(c, "vmax"));
        const std::vector<double> ticks = locator.tickValues(vmin, vmax, space);
        if (!same(ticks, nums(at(c, "ticks")))) {
            failures.add(i, c, "ticks " + show(ticks));
        }
        const Range limits = locator.viewLimits(vmin, vmax);
        if (!same(limits, range(at(c, "view_limits")))) {
            failures.add(i, c, "view_limits " + show(limits));
        }
        ++checked;
    }
    EXPECT_GT(checked, 1000u);
}

TEST(Ticks, MaxNLocatorRejectsNoBins) {
    MaxNLocatorParams params;
    params.nbins = 0;
    EXPECT_THROW(MaxNLocator{params}, std::invalid_argument);
}

TEST(Ticks, FixedLocator) {
    Failures failures("fixed");
    const Array& cases = list(at(data(), "fixed"));
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        std::optional<int> nbins;
        if (has(c, "nbins")) {
            nbins = static_cast<int>(integer(at(c, "nbins")));
        }
        const std::vector<double> locs = nums(at(c, "locs"));
        const std::vector<double> ticks = render::fixedLocatorTicks(locs, nbins);
        if (!same(ticks, nums(at(c, "ticks")))) {
            failures.add(i, c, show(ticks));
        }
    }
}

TEST(Ticks, ScalarFormatter) {
    Failures failures("formatter");
    const Array& cases = list(at(data(), "formatter"));
    ASSERT_GT(cases.size(), 500u);
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        ScalarFormatter formatter(formatterOptions(at(c, "options")));
        const std::vector<double> locs = nums(at(c, "locs"));
        const Range view = range(at(c, "view"));
        if (has(c, "error")) {  // Matplotlib: log10(0) of an empty view interval; MDBridge: exponent 0
            EXPECT_NO_THROW((void)formatter.formatTicks(locs, view));
            continue;
        }
        const std::vector<std::string> labels = formatter.formatTicks(locs, view);
        if (labels != strs(at(c, "labels"))) {
            failures.add(i, c, "labels " + show(labels));
        }
        if (formatter.offsetText() != str(at(c, "offset_text"))) {
            failures.add(i, c, "offset text '" + formatter.offsetText() + "'");
        }
        if (!same(formatter.offset(), num(at(c, "offset"))) ||
            formatter.orderOfMagnitude() != integer(at(c, "order_of_magnitude"))) {
            failures.add(i, c, "offset / order of magnitude");
        }
        if (has(c, "decimals") && formatter.decimals() != integer(at(c, "decimals"))) {
            failures.add(i, c, "decimals " + std::to_string(formatter.decimals()));
        }
    }
}

TEST(Ticks, FormatData) {
    Failures failures("format_data");
    const Array& cases = list(at(data(), "format_data"));
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        ScalarFormatterOptions options;
        options.useMathText = flag(at(c, "mathtext"));
        options.unicodeMinus = flag(at(c, "unicode_minus"));
        const ScalarFormatter formatter(options);
        const std::string text = formatter.formatData(num(at(c, "value")), flag(at(c, "numpy")));
        if (text != str(at(c, "text"))) {
            failures.add(i, c, "'" + text + "'");
        }
    }
}

TEST(Ticks, Axes) {
    Failures failures("axes");
    const Array& cases = list(at(data(), "axes"));
    ASSERT_GT(cases.size(), 100u);
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const std::string mode = str(at(c, "mode"));
        const double labelSize = num(at(c, "label_size"));
        for (const AxisDirection dir : {AxisDirection::X, AxisDirection::Y}) {
            const bool isX = dir == AxisDirection::X;
            const Value& rec = at(c, isX ? "x" : "y");
            Range view;
            if (mode == "limits" || mode == "integer") {
                const Range lim = range(at(c, isX ? "xlim" : "ylim"));
                view = render::setLimits(lim.lo, lim.hi);
            } else {
                const Range dataLim = range(at(c, isX ? "data_x" : "data_y"));
                render::AutoscaleOptions options;
                options.tight = mode == "tight";
                if (mode == "sticky") {
                    options.stickyEdges = nums(at(c, isX ? "sticky_x" : "sticky_y"));
                }
                view = render::autoscale(dataLim.lo, dataLim.hi, MaxNLocator::autoLocator(), options);
            }
            checkAxis(failures, i, c, rec, dir, labelSize, view);
        }
    }
}

TEST(Ticks, SetLimits) {
    Failures failures("set_lim");
    const Array& cases = list(at(data(), "set_lim"));
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const Range r = render::setLimits(num(at(c, "v0")), num(at(c, "v1")));
        if (!same(r, range(at(c, "result")))) {
            failures.add(i, c, show(r));
        }
    }
    EXPECT_THROW((void)render::setLimits(0.0, std::nan("")), std::invalid_argument);
}

TEST(Ticks, Colorbar) {
    Failures failures("colorbar");
    const Array& cases = list(at(data(), "colorbar"));
    ASSERT_GT(cases.size(), 100u);
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const Range limits = render::colorbarLimits(num(at(c, "vmin")), num(at(c, "vmax")));
        if (!same(limits, range(at(c, "limits")))) {
            failures.add(i, c, "limits " + show(limits));
            continue;
        }
        const Value& rec = at(c, "axis");
        checkAxis(failures, i, c, rec, AxisDirection::Y, num(at(c, "font_size")) * 0.75, limits);
        if (!same(nums(at(rec, "ticks")), nums(at(rec, "locs")))) {
            failures.add(i, c, "get_ticks differs from the located ticks");
        }
    }
}

TEST(Ticks, ThinTicks) {
    Failures failures("thin_ticks");
    const Array& cases = list(at(data(), "thin_ticks"));
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const Value& c = cases[i];
        const std::vector<long long> ticks = render::thinTicks(integer(at(c, "n")), integer(at(c, "max_ticks")));
        std::vector<long long> expected;
        for (const Value& x : list(at(c, "ticks"))) {
            expected.push_back(integer(x));
        }
        if (ticks != expected) {
            failures.add(i, c, std::to_string(ticks.size()) + " ticks");
        }
    }
}
