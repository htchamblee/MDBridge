// Axis limits, tick locations and tick labels: a port of Matplotlib 3.11.2's ticker, transforms,
// axis, axes and colorbar code for linear axes (see ticks.hpp), and MDATK's thin_ticks.
//
// Ported from Matplotlib 3.11.2 under Matplotlib's licence (packaging/licenses/Matplotlib-LICENSE.txt):
// Copyright (c) 2012- Matplotlib Development Team; All Rights Reserved.
// Changes made to Matplotlib: translated from Python to C++, linear axes only, rcParams passed as
// options, no locale, no logging of MAXTICKS; the floating-point operations of Python and numpy
// are kept in their order (Python's and numpy's floor division / divmod, numpy's round, C's pow /
// log10) so that the results are Matplotlib's bit for bit.

#include "mdbridge/render/ticks.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <utility>

namespace mdbridge::render {

namespace {

constexpr std::string_view kMinusSign = "\xE2\x88\x92";  // U+2212

// Python's float // and divmod, numpy's floor_divide and divmod for float64 (the same algorithm:
// CPython's float_divmod, numpy's npy_divmod). b == 0 gives numpy's a / b (Python would raise).
double floorDivide(double a, double b, double* modulus = nullptr) {
    if (b == 0.0) {
        if (modulus != nullptr) {
            *modulus = std::fmod(a, b);
        }
        return a / b;
    }
    double mod = std::fmod(a, b);
    double div = (a - mod) / b;
    if (mod != 0.0) {
        if ((b < 0.0) != (mod < 0.0)) {
            mod += b;
            div -= 1.0;
        }
    } else {
        mod = std::copysign(0.0, b);
    }
    double floordiv = 0.0;
    if (div != 0.0) {
        floordiv = std::floor(div);
        if (div - floordiv > 0.5) {
            floordiv += 1.0;
        }
    } else {
        floordiv = std::copysign(0.0, a / b);
    }
    if (modulus != nullptr) {
        *modulus = mod;
    }
    return floordiv;
}

// 10.0 ** k in Python and 10 ** numpy_float in numpy: C's pow(10, k). For whole k the value is
// taken from a table of glibc's results (the correctly rounded power of ten, except for the exact
// ties 10^23 and 10^210, which glibc rounds up), so that other C libraries (MSVC's) give the same.
double cPow10(double k) {
    constexpr int kLow = -340;
    constexpr int kHigh = 308;
    static const std::vector<double> table = [] {
        std::vector<double> out;
        for (int e = kLow; e <= kHigh; ++e) {
            const std::string text = "1e" + std::to_string(e);
            double value = 0.0;
            std::from_chars(text.data(), text.data() + text.size(), value);
            if (e == 23 || e == 210) {
                value = std::nextafter(value, std::numeric_limits<double>::infinity());
            }
            out.push_back(value);
        }
        return out;
    }();
    if (k == std::floor(k) && k >= kLow && k <= kHigh) {
        return table[static_cast<std::size_t>(static_cast<int>(k) - kLow)];
    }
    return std::pow(10.0, k);
}

// The float of Python's int 10 ** k (k >= 0: the exact integer, correctly rounded where it is
// converted to a float; k < 0: Python computes the float power, C's pow).
double pyIntPow10(int k) {
    if (k < 0) {
        return cPow10(static_cast<double>(k));
    }
    const std::string text = "1e" + std::to_string(k);
    double value = 0.0;
    std::from_chars(text.data(), text.data() + text.size(), value);
    return value;
}

// numpy's power_of_ten (numpy/_core/src/multiarray/calculation.c), the factor of np.round.
double numpyPowerOfTen(int n) {
    static constexpr std::array<double, 9> p10 = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8};
    if (n < 9) {
        return p10[static_cast<std::size_t>(n)];
    }
    double ret = 1e9;
    while (n-- > 9) {
        ret *= 10.0;
    }
    return ret;
}

// np.round(x, decimals) of a float64.
double numpyRound(double x, int decimals) {
    if (decimals >= 0) {
        const double f = numpyPowerOfTen(decimals);
        return std::nearbyint(x * f) / f;
    }
    const double f = numpyPowerOfTen(-decimals);
    return std::nearbyint(x / f) * f;
}

// '%1.<precision>f' % x (correctly rounded, as Python's float formatting).
std::string formatFixed(double x, int precision) {
    if (std::isnan(x)) {
        return "nan";
    }
    if (std::isinf(x)) {
        return x < 0 ? "-inf" : "inf";
    }
    std::array<char, 400> buffer{};
    auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), x, std::chars_format::fixed, precision);
    if (result.ec != std::errc()) {  // more than 400 characters: precision beyond any tick label
        return std::to_string(x);
    }
    return std::string(buffer.data(), result.ptr);
}

// Python's round(x, ndigits) of a float: x rounded correctly to ndigits decimals.
double pythonRound(double x, int ndigits) {
    if (!std::isfinite(x)) {
        return x;
    }
    const std::string text = formatFixed(x, ndigits);
    double value = 0.0;
    std::from_chars(text.data(), text.data() + text.size(), value);
    return std::copysign(value, x);
}

// '%1.10g' % x
std::string formatGeneral10(double x) {
    std::array<char, 64> buffer{};
    auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), x, std::chars_format::general, 10);
    return std::string(buffer.data(), result.ptr);
}

Range sortedRange(Range view) {
    if (view.hi < view.lo) {
        std::swap(view.lo, view.hi);
    }
    return view;
}

std::vector<double> visible(std::span<const double> locs, Range view) {
    const Range v = sortedRange(view);
    std::vector<double> out;
    for (const double x : locs) {
        if (v.lo <= x && x <= v.hi) {
            out.push_back(x);
        }
    }
    return out;
}

std::vector<double> validateSteps(const std::vector<double>& steps) {
    static const char* message = "steps argument must be an increasing sequence of numbers between 1 and 10 inclusive";
    if (steps.empty()) {
        throw std::invalid_argument(message);
    }
    for (std::size_t i = 1; i < steps.size(); ++i) {
        if (!(steps[i] - steps[i - 1] > 0)) {
            throw std::invalid_argument(message);
        }
    }
    if (steps.back() > 10 || steps.front() < 1) {
        throw std::invalid_argument(message);
    }
    std::vector<double> out;
    if (steps.front() != 1) {
        out.push_back(1.0);
    }
    out.insert(out.end(), steps.begin(), steps.end());
    if (out.back() != 10) {
        out.push_back(10.0);
    }
    return out;
}

// _Edge_integer
class EdgeInteger {
public:
    EdgeInteger(double step, double offset) : m_step(step), m_offset(std::abs(offset)) {}

    bool closeTo(double ms, double edge) const {
        double tol = 1e-10;
        if (m_offset > 0) {
            const double digits = std::log10(m_offset / m_step);
            tol = std::max(1e-10, cPow10(digits - 12));
            tol = std::min(0.4999, tol);
        }
        return std::abs(ms - edge) < tol;
    }

    double le(double x) const {
        double m = 0.0;
        const double d = floorDivide(x, m_step, &m);
        return closeTo(m / m_step, 1) ? d + 1 : d;
    }

    double ge(double x) const {
        double m = 0.0;
        const double d = floorDivide(x, m_step, &m);
        return closeTo(m / m_step, 0) ? d : d + 1;
    }

private:
    double m_step;
    double m_offset;
};

}  // namespace

Range nonsingular(double vmin, double vmax, double expander, double tiny, bool increasing) {
    if (!std::isfinite(vmin) || !std::isfinite(vmax)) {
        return {-expander, expander};
    }
    bool swapped = false;
    if (vmax < vmin) {
        std::swap(vmin, vmax);
        swapped = true;
    }
    const double maxabsvalue = std::max(std::abs(vmin), std::abs(vmax));
    if (maxabsvalue < (1e6 / tiny) * std::numeric_limits<double>::min()) {
        vmin = -expander;
        vmax = expander;
    } else if (vmax - vmin <= maxabsvalue * tiny) {
        if (vmax == 0 && vmin == 0) {
            vmin = -expander;
            vmax = expander;
        } else {
            vmin -= expander * std::abs(vmin);
            vmax += expander * std::abs(vmax);
        }
    }
    if (swapped && !increasing) {
        std::swap(vmin, vmax);
    }
    return {vmin, vmax};
}

double fontSizeInPoints(double size) { return size < 1.0 ? 1.0 : size; }

int tickSpace(AxisDirection direction, double lengthPoints, double labelSizePoints) {
    const double size = fontSizeInPoints(labelSizePoints) * (direction == AxisDirection::X ? 3 : 2);
    if (!(size > 0)) {
        return 2147483647;
    }
    const double n = std::floor(lengthPoints / size);
    // int(np.floor(...)) of an axis length: far inside int's range for any real figure
    if (n >= 2147483647.0) {
        return 2147483647;
    }
    if (n <= -2147483647.0) {
        return -2147483647;
    }
    return static_cast<int>(n);
}

// ---------------------------------------------------------------------------------------------
// MaxNLocator

MaxNLocator::MaxNLocator(MaxNLocatorParams params) : m_params(std::move(params)) {
    m_nbins = m_params.nbins;
    if (m_nbins && *m_nbins < 1) {
        throw std::invalid_argument("nbins must be at least 1 or 'auto'");
    }
    m_params.minNTicks = std::max(1, m_params.minNTicks);
    std::vector<double> steps = m_params.steps.empty()
                                    ? std::vector<double>{1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10}
                                    : validateSteps(m_params.steps);
    m_params.steps = steps;
    // _staircase: 0.1 x steps[:-1], steps, 10 x steps[1]
    for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
        m_extendedSteps.push_back(0.1 * steps[i]);
    }
    m_extendedSteps.insert(m_extendedSteps.end(), steps.begin(), steps.end());
    m_extendedSteps.push_back(10 * steps[1]);
}

MaxNLocator MaxNLocator::autoLocator() {
    MaxNLocatorParams p;
    p.nbins = std::nullopt;
    p.steps = {1, 2, 2.5, 5, 10};
    return MaxNLocator(std::move(p));
}

MaxNLocator MaxNLocator::integerLocator() {
    MaxNLocatorParams p;
    p.integer = true;
    return MaxNLocator(std::move(p));
}

int MaxNLocator::bins(std::optional<int> tickSpace) const {
    if (m_nbins) {
        return *m_nbins;
    }
    if (!tickSpace) {
        return 9;
    }
    return std::clamp(*tickSpace, std::max(1, m_params.minNTicks - 1), 9);  // np.clip
}

std::vector<double> MaxNLocator::rawTicks(double vmin, double vmax, std::optional<int> tickSpace) const {
    const int nbins = bins(tickSpace);
    const double n = static_cast<double>(nbins);

    // scale_range(vmin, vmax, nbins)
    const double dv = std::abs(vmax - vmin);
    const double meanv = (vmax + vmin) / 2;
    double offset = 0.0;
    if (!(std::abs(meanv) / dv < 100)) {
        offset = std::copysign(cPow10(floorDivide(std::log10(std::abs(meanv)), 1.0)), meanv);
    }
    const double scale = cPow10(floorDivide(std::log10(dv / n), 1.0));

    const double vmin0 = vmin - offset;
    const double vmax0 = vmax - offset;
    std::vector<double> steps;
    steps.reserve(m_extendedSteps.size());
    for (const double s : m_extendedSteps) {
        const double step = s * scale;
        if (m_params.integer && !(step < 1 || std::abs(step - std::nearbyint(step)) < 0.001)) {
            continue;  // for steps > 1, keep only integer values
        }
        steps.push_back(step);
    }
    const double rawStep = (vmax0 - vmin0) / n;

    // index of the smallest step >= the raw step (else the largest)
    std::size_t istep = steps.empty() ? 0 : steps.size() - 1;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (steps[i] >= rawStep) {
            istep = i;
            break;
        }
    }

    std::vector<double> ticks;
    for (std::size_t k = std::min(istep + 1, steps.size()); k-- > 0;) {
        double step = steps[k];
        if (m_params.integer &&
            std::floor(vmax0) - std::ceil(vmin0) >= static_cast<double>(m_params.minNTicks - 1)) {
            step = std::max(1.0, step);
        }
        const double bestVmin = floorDivide(vmin0, step) * step;
        if (!(step > 0)) {
            throw std::invalid_argument("'step' must be positive");
        }
        const EdgeInteger edge(step, offset);
        const double low = edge.le(vmin0 - bestVmin);
        const double high = edge.ge(vmax0 - bestVmin);
        // np.arange(low, high + 1) * step + best_vmin
        ticks.clear();
        const double count = std::ceil(high + 1 - low);
        for (double i = 0; i < count; i += 1) {
            ticks.push_back((low + i) * step + bestVmin);
        }
        std::ptrdiff_t nticks = 0;
        for (const double t : ticks) {
            if (t <= vmax0 && t >= vmin0) {
                ++nticks;
            }
        }
        if (nticks >= m_params.minNTicks) {
            break;
        }
    }
    for (double& t : ticks) {
        t += offset;
    }
    return ticks;
}

std::vector<double> MaxNLocator::tickValues(double vmin, double vmax, std::optional<int> tickSpace) const {
    if (m_params.symmetric) {
        vmax = std::max(std::abs(vmin), std::abs(vmax));
        vmin = -vmax;
    }
    const Range r = nonsingular(vmin, vmax, 1e-13, 1e-14);
    std::vector<double> locs = rawTicks(r.lo, r.hi, tickSpace);
    auto dropFront = [&locs] {
        if (!locs.empty()) {
            locs.erase(locs.begin());
        }
    };
    auto dropBack = [&locs] {
        if (!locs.empty()) {
            locs.pop_back();
        }
    };
    switch (m_params.prune) {
        case Prune::Lower: dropFront(); break;
        case Prune::Upper: dropBack(); break;
        case Prune::Both: dropFront(); dropBack(); break;
        case Prune::None: break;
    }
    return locs;
}

Range MaxNLocator::viewLimits(double dmin, double dmax) const {
    if (m_params.symmetric) {
        dmax = std::max(std::abs(dmin), std::abs(dmax));
        dmin = -dmax;
    }
    return nonsingular(dmin, dmax, 1e-12, 1e-13);
}

Range MaxNLocator::nonsingularLimits(double v0, double v1) { return nonsingular(v0, v1, 0.05); }

std::vector<double> fixedLocatorTicks(std::span<const double> locs, std::optional<int> nbins) {
    if (!nbins || locs.empty()) {
        return {locs.begin(), locs.end()};
    }
    nbins = std::max(*nbins, 2);  // FixedLocator.__init__
    const auto n = static_cast<std::ptrdiff_t>(locs.size());
    const auto step = static_cast<std::ptrdiff_t>(
        std::max(std::ceil(static_cast<double>(n) / static_cast<double>(*nbins)), 1.0));
    auto every = [&](std::ptrdiff_t start) {
        std::vector<double> out;
        for (std::ptrdiff_t i = start; i < n; i += step) {
            out.push_back(locs[static_cast<std::size_t>(i)]);
        }
        return out;
    };
    auto minAbs = [](const std::vector<double>& v) {
        double m = std::numeric_limits<double>::infinity();
        for (const double x : v) {
            m = std::min(m, std::abs(x));
        }
        return m;
    };
    std::vector<double> ticks = every(0);
    for (std::ptrdiff_t i = 1; i < step; ++i) {
        std::vector<double> other = every(i);
        if (!other.empty() && minAbs(other) < minAbs(ticks)) {
            ticks = std::move(other);
        }
    }
    return ticks;
}

// ---------------------------------------------------------------------------------------------
// ScalarFormatter

ScalarFormatter::ScalarFormatter(ScalarFormatterOptions options) : m_options(std::move(options)) {
    if (m_options.fixedOffset) {
        m_options.useOffset = false;
        m_offset = *m_options.fixedOffset;
        m_offsetFromNumpy = false;
    }
}

std::string ScalarFormatter::fixMinus(std::string text) const {
    if (!m_options.unicodeMinus) {
        return text;
    }
    std::string out;
    out.reserve(text.size() + 8);
    for (const char c : text) {
        if (c == '-') {
            out += kMinusSign;
        } else {
            out += c;
        }
    }
    return out;
}

void ScalarFormatter::setLocs(std::span<const double> locs, Range view) {
    m_locs.assign(locs.begin(), locs.end());
    if (m_locs.empty()) {
        return;
    }
    if (m_options.useOffset) {
        computeOffset(view);
    }
    setOrderOfMagnitude(view);
    setFormat(view);
}

void ScalarFormatter::computeOffset(Range view) {
    m_offsetFromNumpy = true;
    const std::vector<double> locs = visible(m_locs, view);
    if (locs.empty()) {
        m_offset = 0;
        return;
    }
    const double lmin = *std::min_element(locs.begin(), locs.end());
    const double lmax = *std::max_element(locs.begin(), locs.end());
    if (lmin == lmax || (lmin <= 0 && 0 <= lmax)) {
        m_offset = 0;
        return;
    }
    double absMin = std::abs(lmin);
    double absMax = std::abs(lmax);
    if (absMax < absMin) {
        std::swap(absMin, absMax);
    }
    const double sign = std::copysign(1.0, lmin);
    const double oomMax = std::ceil(std::log10(absMax));
    // the oom loops end long before 10 ** oom underflows; the bound only guards against a hang
    constexpr double kLowest = -400;
    double oom = oomMax;
    while (oom > kLowest && floorDivide(absMin, cPow10(oom)) == floorDivide(absMax, cPow10(oom))) {
        oom -= 1;
    }
    oom += 1;
    if ((absMax - absMin) / cPow10(oom) <= 1e-2) {
        oom = oomMax;
        while (oom > kLowest && !(floorDivide(absMax, cPow10(oom)) - floorDivide(absMin, cPow10(oom)) > 1)) {
            oom -= 1;
        }
        oom += 1;
    }
    const int n = m_options.offsetThreshold - 1;
    const double q = floorDivide(absMax, cPow10(oom));
    m_offset = q >= pyIntPow10(n) ? sign * q * cPow10(oom) : 0.0;
}

void ScalarFormatter::setOrderOfMagnitude(Range view) {
    if (!m_options.scientific) {
        m_orderOfMagnitude = 0;
        return;
    }
    if (m_options.powerLimitLow == m_options.powerLimitHigh && m_options.powerLimitLow != 0) {
        m_orderOfMagnitude = m_options.powerLimitLow;
        return;
    }
    const Range v = sortedRange(view);
    const std::vector<double> locs = visible(m_locs, v);
    if (locs.empty()) {
        m_orderOfMagnitude = 0;
        return;
    }
    int oom = 0;
    if (m_offset != 0) {
        // Matplotlib raises a math domain error for an empty view interval; 0 here
        const double span = v.hi - v.lo;
        oom = span > 0 ? static_cast<int>(std::floor(std::log10(span))) : 0;
    } else {
        double val = 0;
        for (const double x : locs) {
            val = std::max(val, std::abs(x));
        }
        oom = val == 0 ? 0 : static_cast<int>(std::floor(std::log10(val)));
    }
    if (oom <= m_options.powerLimitLow || oom >= m_options.powerLimitHigh) {
        m_orderOfMagnitude = oom;
    } else {
        m_orderOfMagnitude = 0;
    }
}

void ScalarFormatter::setFormat(Range view) {
    std::vector<double> raw = m_locs;
    const bool few = m_locs.size() < 2;
    if (few) {
        raw.push_back(view.lo);
        raw.push_back(view.hi);
    }
    const double scale = cPow10(static_cast<double>(m_orderOfMagnitude));
    std::vector<double> locs;
    locs.reserve(raw.size());
    for (const double x : raw) {
        locs.push_back((x - m_offset) / scale);
    }
    double lo = locs.front();
    double hi = locs.front();
    for (const double x : locs) {
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    double locRange = hi - lo;  // np.ptp
    if (locRange == 0) {
        locRange = 0;
        for (const double x : locs) {
            locRange = std::max(locRange, std::abs(x));
        }
    }
    if (locRange == 0) {
        locRange = 1;
    }
    if (few) {
        locs.resize(locs.size() - 2);
    }
    const int locRangeOom = static_cast<int>(std::floor(std::log10(locRange)));
    int sigfigs = std::max(0, 3 - locRangeOom);
    const double thresh = 1e-3 * pyIntPow10(locRangeOom);
    while (sigfigs >= 0) {
        double maxDev = -std::numeric_limits<double>::infinity();
        for (const double x : locs) {
            maxDev = std::max(maxDev, std::abs(x - numpyRound(x, sigfigs)));
        }
        if (maxDev < thresh) {
            --sigfigs;
        } else {
            break;
        }
    }
    m_decimals = sigfigs + 1;
    m_hasFormat = true;
}

std::string ScalarFormatter::operator()(double value) const {
    if (m_locs.empty()) {
        return "";
    }
    double xp = (value - m_offset) / cPow10(static_cast<double>(m_orderOfMagnitude));
    if (std::abs(xp) < 1e-8) {
        xp = 0;  // Python's int 0: '0.00', never '-0.00'
    }
    std::string text = formatFixed(xp, m_hasFormat ? m_decimals : 0);
    if (m_options.useMathText) {
        text = "$\\mathdefault{" + text + "}$";
    }
    return fixMinus(std::move(text));
}

std::vector<std::string> ScalarFormatter::formatTicks(std::span<const double> locs, Range view) {
    setLocs(locs, view);
    std::vector<std::string> out;
    out.reserve(locs.size());
    for (const double x : locs) {
        out.push_back((*this)(x));
    }
    return out;
}

std::string ScalarFormatter::formatData(double value, bool numpyValue) const {
    const int e = static_cast<int>(std::floor(std::log10(std::abs(value))));
    const double divided = value / pyIntPow10(e);
    const double s = numpyValue ? numpyRound(divided, 10) : pythonRound(divided, 10);
    std::string significand;
    if (std::fmod(s, 1.0) == 0) {
        significand = formatFixed(s, 0);  // '%d' of an integral float
        if (significand == "-0") {
            significand = "0";
        }
    } else {
        significand = formatGeneral10(s);
    }
    significand = fixMinus(std::move(significand));
    if (e == 0) {
        return significand;
    }
    std::string exponent = fixMinus(std::to_string(e));
    if (m_options.useMathText) {
        exponent = "10^{" + exponent + "}";
        return s == 1 ? exponent : significand + " \\times " + exponent;
    }
    return significand + "e" + exponent;
}

std::string ScalarFormatter::offsetText() const {
    if (m_locs.empty()) {
        return "";
    }
    if (m_orderOfMagnitude == 0 && m_offset == 0) {
        return "";
    }
    std::string offsetStr;
    std::string sciNotStr;
    if (m_offset != 0) {
        offsetStr = formatData(m_offset, m_offsetFromNumpy);
        if (m_offset > 0) {
            offsetStr = "+" + offsetStr;
        }
    }
    if (m_orderOfMagnitude != 0) {
        if (m_options.useMathText) {
            // format_data(10 ** oom): a Python int for oom >= 0, else a float
            sciNotStr = formatData(pyIntPow10(m_orderOfMagnitude), false);
        } else {
            sciNotStr = "1e" + std::to_string(m_orderOfMagnitude);
        }
    }
    std::string s;
    if (m_options.useMathText) {
        if (!sciNotStr.empty()) {
            sciNotStr = "\\times\\mathdefault{" + sciNotStr + "}";
        }
        s = "$" + sciNotStr + "\\mathdefault{" + offsetStr + "}$";
    } else {
        s = sciNotStr + offsetStr;
    }
    return fixMinus(std::move(s));
}

// ---------------------------------------------------------------------------------------------
// Limits

Range setLimits(double v0, double v1) {
    if (!std::isfinite(v0) || !std::isfinite(v1)) {
        throw std::invalid_argument("Axis limits cannot be NaN or Inf");
    }
    const bool reverse = v0 > v1;
    Range r = MaxNLocator::nonsingularLimits(v0, v1);
    if (r.hi < r.lo) {
        std::swap(r.lo, r.hi);
    }
    if (reverse) {
        std::swap(r.lo, r.hi);
    }
    return r;
}

Range autoscale(double dataMin, double dataMax, const MaxNLocator& locator, const AutoscaleOptions& options) {
    const double inf = std::numeric_limits<double>::infinity();
    double x0 = -inf;
    double x1 = inf;
    const bool f0 = std::isfinite(dataMin);
    const bool f1 = std::isfinite(dataMax);
    if (f0 || f1) {
        x0 = f0 && f1 ? std::min(dataMin, dataMax) : (f0 ? dataMin : dataMax);
        x1 = f0 && f1 ? std::max(dataMin, dataMax) : (f0 ? dataMin : dataMax);
    }
    Range r = MaxNLocator::nonsingularLimits(x0, x1);
    x0 = r.lo;
    x1 = r.hi;

    std::optional<double> x0bound;
    std::optional<double> x1bound;
    if (options.margin != 0 && !options.stickyEdges.empty()) {
        std::vector<double> stickies = options.stickyEdges;
        std::sort(stickies.begin(), stickies.end());
        const double tol = 1e-5 * std::abs(x1 - x0);
        // index of the largest element < x0 + tol, of the smallest element > x1 - tol
        const auto i0 = std::lower_bound(stickies.begin(), stickies.end(), x0 + tol) - stickies.begin() - 1;
        if (i0 != -1) {
            x0bound = stickies[static_cast<std::size_t>(i0)];
        }
        const auto i1 = std::lower_bound(stickies.begin(), stickies.end(), x1 - tol) - stickies.begin();
        if (i1 != static_cast<std::ptrdiff_t>(stickies.size())) {
            x1bound = stickies[static_cast<std::size_t>(i1)];
        }
    }
    double delta = (x1 - x0) * options.margin;
    if (!std::isfinite(delta)) {
        delta = 0;
    }
    x0 = x0 - delta;
    x1 = x1 + delta;
    if (x0bound) {
        x0 = std::max(x0, *x0bound);
    }
    if (x1bound) {
        x1 = std::min(x1, *x1bound);
    }
    if (!options.tight) {
        r = locator.viewLimits(x0, x1);
        x0 = r.lo;
        x1 = r.hi;
    }
    // set_xbound: sorted, reversed on an inverted axis; then set_xlim
    if (x1 < x0) {
        std::swap(x0, x1);
    }
    if (options.inverted) {
        std::swap(x0, x1);
    }
    return setLimits(x0, x1);
}

std::vector<std::size_t> ticksToDraw(std::span<const double> locs, Range view) {
    const Range v = sortedRange(view);
    const double rtol = (v.hi - v.lo) * 1e-10;
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < locs.size(); ++i) {
        if (v.lo - rtol <= locs[i] && locs[i] <= v.hi + rtol) {
            out.push_back(i);
        }
    }
    return out;
}

Range colorbarLimits(double vmin, double vmax, bool inverted) {
    const Range norm = nonsingular(vmin, vmax, 0.1);
    // Normalize.inverse of the boundaries 0 and 1 (np.linspace keeps both ends exact)
    const double lower = norm.lo + 0.0 * (norm.hi - norm.lo);
    const double upper = norm.lo + 1.0 * (norm.hi - norm.lo);
    return inverted ? setLimits(upper, lower) : setLimits(lower, upper);
}

std::vector<long long> thinTicks(long long n, long long maxTicks) {
    maxTicks = std::max<long long>(1, maxTicks);
    const long long step = n > 0 ? std::max<long long>(1, (n + maxTicks - 1) / maxTicks) : 1;
    std::vector<long long> out;
    for (long long i = 0; i < n; i += step) {
        out.push_back(i);
    }
    return out;
}

}  // namespace mdbridge::render
