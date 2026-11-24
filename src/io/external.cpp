// numpyInterp is a port of numpy's arr_interp / binary_search_with_guess
// (numpy/_core/src/multiarray/compiled_base.c), used under numpy's licence:
//
// Copyright (c) 2005-2025, NumPy Developers.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice, this list of
//       conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright notice, this list of
//       conditions and the following disclaimer in the documentation and/or other materials
//       provided with the distribution.
//
//     * Neither the name of the NumPy Developers nor the names of any contributors may be used to
//       endorse or promote products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "mdbridge/io/external.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"
#include "text_file.hpp"

namespace mdbridge::io {

namespace {

using core::ValueError;

std::string utf8(std::u32string_view text) { return core::encodeUtf8(text); }

std::u32string strip(std::u32string_view text) {
    std::size_t first = 0;
    std::size_t last = text.size();
    while (first < last && core::isPythonSpace(text[first])) {
        ++first;
    }
    while (last > first && core::isPythonSpace(text[last - 1])) {
        --last;
    }
    return std::u32string(text.substr(first, last - first));
}

// str.split(sep)
std::vector<std::u32string> splitOn(std::u32string_view text, std::u32string_view sep) {
    std::vector<std::u32string> out;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t at = text.find(sep, begin);
        if (at == std::u32string_view::npos) {
            out.emplace_back(text.substr(begin));
            return out;
        }
        out.emplace_back(text.substr(begin, at - begin));
        begin = at + sep.size();
    }
}

// str.split()
std::vector<std::u32string> splitWhitespace(std::u32string_view text) {
    std::vector<std::u32string> out;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && core::isPythonSpace(text[i])) {
            ++i;
        }
        const std::size_t start = i;
        while (i < text.size() && !core::isPythonSpace(text[i])) {
            ++i;
        }
        if (i > start) {
            out.emplace_back(text.substr(start, i - start));
        }
    }
    return out;
}

// _split: the explicit delimiter, else ",", else tab, else whitespace
std::vector<std::string> splitFields(std::u32string_view line, std::u32string_view delimiter) {
    std::vector<std::u32string> parts;
    bool stripParts = true;
    if (!delimiter.empty()) {
        parts = splitOn(line, delimiter);
    } else if (line.find(U',') != std::u32string_view::npos) {
        parts = splitOn(line, U",");
    } else if (line.find(U'\t') != std::u32string_view::npos) {
        parts = splitOn(line, U"\t");
    } else {
        parts = splitWhitespace(line);
        stripParts = false;
    }
    std::vector<std::string> out;
    out.reserve(parts.size());
    for (const std::u32string& part : parts) {
        out.push_back(utf8(stripParts ? strip(part) : part));
    }
    return out;
}

// numpy's binary_search_with_guess: the index i with arr[i] <= key < arr[i + 1], -1 below arr[0]
// and len above arr[len - 1], starting from a guess (the previous result) since the keys usually
// come in order.
constexpr std::ptrdiff_t kLikelyInCacheSize = 8;

std::ptrdiff_t binarySearchWithGuess(double key, const double* arr, std::ptrdiff_t len,
                                     std::ptrdiff_t guess) {
    std::ptrdiff_t imin = 0;
    std::ptrdiff_t imax = len;

    if (key > arr[len - 1]) {
        return len;
    }
    if (key < arr[0]) {
        return -1;
    }

    // linear search for short arrays; from above key >= arr[0]
    if (len <= 4) {
        std::ptrdiff_t i = 1;
        for (; i < len && key >= arr[i]; ++i) {
        }
        return i - 1;
    }

    if (guess > len - 3) {
        guess = len - 3;
    }
    if (guess < 1) {
        guess = 1;
    }

    // check the most likely values: guess - 1, guess, guess + 1
    if (key < arr[guess]) {
        if (key < arr[guess - 1]) {
            imax = guess - 1;
            // last attempt to restrict the search to items in cache
            if (guess > kLikelyInCacheSize && key >= arr[guess - kLikelyInCacheSize]) {
                imin = guess - kLikelyInCacheSize;
            }
        } else {
            return guess - 1;
        }
    } else {
        if (key < arr[guess + 1]) {
            return guess;
        }
        if (key < arr[guess + 2]) {
            return guess + 1;
        }
        imin = guess + 2;
        if (guess < len - kLikelyInCacheSize - 1 && key < arr[guess + kLikelyInCacheSize]) {
            imax = guess + kLikelyInCacheSize;
        }
    }

    // bisection
    while (imin < imax) {
        const std::ptrdiff_t imid = imin + ((imax - imin) >> 1);
        if (key >= arr[imid]) {
            imin = imid + 1;
        } else {
            imax = imid;
        }
    }
    return imin - 1;
}

}  // namespace

DataTable readTable(const std::filesystem::path& path, std::string_view delimiter) {
    const std::u32string delim = core::decodeUtf8(delimiter);
    std::vector<std::string> header;
    std::vector<std::vector<double>> rows;
    std::vector<std::string> pendingComment;
    for (const std::u32string& raw : detail::readTextLines(path)) {
        const std::u32string line = strip(raw);
        if (line.empty()) {
            continue;
        }
        if (std::u32string_view(U"#@;%").find(line[0]) != std::u32string_view::npos) {
            if (rows.empty()) {
                const std::size_t start = line.find_first_not_of(U"#@;% ");
                pendingComment =
                    splitFields(start == std::u32string::npos ? U"" : line.substr(start), delim);
            }
            continue;
        }
        const std::vector<std::string> tokens = splitFields(line, delim);
        if (rows.empty() && header.empty() &&
            !std::all_of(tokens.begin(), tokens.end(), [](const std::string& t) {
                return core::tryPythonFloat(t).has_value();
            })) {
            header = tokens;
            continue;
        }
        std::vector<double> row;
        for (const std::string& token : tokens) {
            const std::optional<double> value = core::tryPythonFloat(token);
            if (!value) {
                throw ValueError(detail::pathText(path) + ": non-numeric data line " +
                                 core::pythonStringRepr(utf8(line)));
            }
            row.push_back(*value);
        }
        rows.push_back(std::move(row));
    }
    if (rows.empty()) {
        throw ValueError(detail::pathText(path) + ": no numeric data found");
    }
    std::size_t width = rows.front().size();
    for (const auto& row : rows) {
        width = std::min(width, row.size());
    }
    for (auto& row : rows) {
        row.resize(width);
    }
    if (header.empty() && pendingComment.size() == width) {
        header = pendingComment;
    }
    if (header.size() != width) {
        header.clear();
        for (std::size_t k = 0; k < width; ++k) {
            header.push_back("col" + std::to_string(k));
        }
    }
    return {std::move(header), std::move(rows)};
}

std::vector<double> numpyInterp(std::span<const double> target, std::span<const double> x,
                                std::span<const double> y, double left, double right) {
    const auto lenxp = static_cast<std::ptrdiff_t>(x.size());
    if (lenxp == 0 && !target.empty()) {
        throw ValueError("array of sample points is empty");
    }
    if (y.size() != x.size()) {
        throw ValueError("fp and xp are not of the same length.");
    }
    std::vector<double> out(target.size());
    if (target.empty()) {
        return out;
    }
    const double* dx = x.data();
    const double* dy = y.data();
    if (lenxp == 1) {
        // binary_search_with_guess needs at least a 3 item long array
        for (std::size_t i = 0; i < target.size(); ++i) {
            const double xv = target[i];
            out[i] = xv < dx[0] ? left : (xv > dx[0] ? right : dy[0]);
        }
        return out;
    }
    std::ptrdiff_t j = 0;
    for (std::size_t i = 0; i < target.size(); ++i) {
        const double xv = target[i];
        if (std::isnan(xv)) {
            out[i] = xv;
            continue;
        }
        j = binarySearchWithGuess(xv, dx, lenxp, j);
        if (j == -1) {
            out[i] = left;
        } else if (j == lenxp) {
            out[i] = right;
        } else if (j == lenxp - 1) {
            out[i] = dy[j];
        } else if (dx[j] == xv) {
            out[i] = dy[j];  // avoids a non-finite interpolation
        } else {
            const double slope = (dy[j + 1] - dy[j]) / (dx[j + 1] - dx[j]);
            // if one direction gives NaN, try the other
            const double fromLeft = slope * (xv - dx[j]);
            out[i] = fromLeft + dy[j];
            if (std::isnan(out[i])) {
                const double fromRight = slope * (xv - dx[j + 1]);
                out[i] = fromRight + dy[j + 1];
                if (std::isnan(out[i]) && dy[j] == dy[j + 1]) {
                    out[i] = dy[j];
                }
            }
        }
    }
    return out;
}

int columnIndex(const std::vector<std::string>& header, std::string_view spec, int fallback) {
    if (spec.empty()) {
        return fallback;
    }
    if (const auto it = std::find(header.begin(), header.end(), spec); it != header.end()) {
        return static_cast<int>(it - header.begin());
    }
    const std::u32string code = core::decodeUtf8(spec);
    if (std::all_of(code.begin(), code.end(), core::isPythonDigit)) {
        try {
            const std::int64_t k = core::pythonInt(spec);
            if (k < static_cast<std::int64_t>(header.size())) {
                return static_cast<int>(k);
            }
        } catch (const std::overflow_error&) {
            // beyond 64 bits: no such column
        }
    }
    std::string available;
    for (std::size_t k = 0; k < header.size(); ++k) {
        available += (k ? ", " : "") + header[k];
    }
    throw ValueError("column " + core::pythonStringRepr(std::string(spec)) +
                     " not found (available: " + available + ")");
}

std::string sanitizeKey(std::string_view text) {
    // The allowed characters are ASCII, so every byte of another character is outside them and a
    // run of such bytes becomes one "_", as the run of characters does in Python.
    std::string out;
    bool inRun = false;
    for (const char c : text) {
        const bool allowed = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                             (c >= '0' && c <= '9') || c == '_' || c == '.' || c == '+' || c == '-';
        if (allowed) {
            out += c;
        } else if (!inRun) {
            out += '_';
        }
        inRun = !allowed;
    }
    const std::size_t first = out.find_first_not_of('_');
    if (first == std::string::npos) {
        return "value";
    }
    return out.substr(first, out.find_last_not_of('_') - first + 1);
}

core::AnalysisResult loadExternal(const std::vector<core::ExternalSeries>& specs,
                                  const core::Dataset& dataset) {
    core::AnalysisResult result;
    result.name = "external";
    for (std::size_t n = 0; n < specs.size(); ++n) {
        const core::ExternalSeries& spec = specs[n];
        if (spec.path.empty()) {
            continue;
        }
        const DataTable table = readTable(detail::pathFromUtf8(spec.path), spec.delimiter);
        const std::vector<std::string>& header = table.header;
        const int width = static_cast<int>(header.size());
        const auto column = [&](int index) {
            const auto c = static_cast<std::size_t>(index < 0 ? index + width : index);
            std::vector<double> values;
            values.reserve(table.rows.size());
            for (const auto& row : table.rows) {
                values.push_back(row.at(c));
            }
            return values;
        };
        const int xi = columnIndex(header, spec.x_column, 0);
        std::vector<double> x = column(xi);
        const std::string unit = spec.x_unit.empty() ? "ps" : spec.x_unit;
        std::vector<double> target;
        if (unit == "frame" || unit == "frame1") {
            for (const std::int64_t f : dataset.frames()) {
                target.push_back(static_cast<double>(f) + (unit == "frame1" ? 1.0 : 0.0));
            }
        } else {
            const double factor = core::timeUnitFactor(unit);
            for (double& v : x) {
                v *= factor;
            }
            target = dataset.timePs();
        }
        // np.argsort(x, kind="stable"): NaN last
        std::vector<std::size_t> order(x.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return x[a] < x[b] || (!std::isnan(x[a]) && std::isnan(x[b]));
        });
        std::vector<double> xs;
        for (const std::size_t k : order) {
            xs.push_back(x[k]);
        }
        std::vector<std::string> columns = spec.columns;
        if (columns.empty()) {
            for (int k = 0; k < width; ++k) {
                if (k != xi) {
                    columns.push_back(header[static_cast<std::size_t>(k)]);
                }
            }
        }
        const std::string label =
            !spec.label.empty()
                ? spec.label
                : core::posixpath::splitext(core::pathBasename(spec.path)).first;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (const std::string& col : columns) {
            const int ci = columnIndex(header, col, -1);
            const std::vector<double> data = column(ci);
            std::vector<double> y;
            for (const std::size_t k : order) {
                y.push_back(data[k]);
            }
            const std::string& name = header[static_cast<std::size_t>(ci < 0 ? ci + width : ci)];
            const std::string array = "s" + std::to_string(n) + "_" + std::to_string(ci);
            result.setArray(array, core::NumArray(numpyInterp(target, xs, y, nan, nan)));
            result.addSeries(sanitizeKey(label + "." + name), array, label + ": " + name, "");
        }
    }
    return result;
}

}  // namespace mdbridge::io
