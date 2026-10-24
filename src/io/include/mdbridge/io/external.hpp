#pragma once
// External time series (mdatk/io/external.py): CSV or whitespace text tables such as cpptraj .dat
// files, interpolated onto the analysed frames as extra series.

#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/project.hpp"

namespace mdbridge::io {

struct DataTable {
    std::vector<std::string> header;
    std::vector<std::vector<double>> rows;  // all rows cut to the shortest
};

// read_table(path, delimiter): comment lines (# @ ; %), the header line or the last comment before
// the data, delimiters detected per line. Throws core::ValueError ("<path>: no numeric data
// found", "<path>: non-numeric data line '...'"), std::filesystem::filesystem_error.
DataTable readTable(const std::filesystem::path& path, std::string_view delimiter = {});

// numpy.interp(target, x, y, left, right), numpy's algorithm exactly.
std::vector<double> numpyInterp(std::span<const double> target, std::span<const double> x, std::span<const double> y,
                                double left, double right);

// The column `spec` names (header name first, then a 0-based index); `fallback` for "".
// Throws core::ValueError("column 'x' not found (available: a, b)").
int columnIndex(const std::vector<std::string>& header, std::string_view spec, int fallback);

// Output-key cleaner: runs outside [A-Za-z0-9_.+-] become "_", "_" stripped, "value" if empty.
std::string sanitizeKey(std::string_view text);

// load_external(specs, dataset): result "external" with one series per column.
core::AnalysisResult loadExternal(const std::vector<core::ExternalSeries>& specs, const core::Dataset& dataset);

}  // namespace mdbridge::io
