#pragma once
// Analysis results as MDATK hands them to panels, charts and exports (mdatk/data.py): one result
// per analysis with its arrays and the description of its outputs, and the dataset of all results
// over the analysed frames with their time axis.
//
// Lookups that fail throw std::out_of_range with MDATK's KeyError messages, e.g. "analysis 'x' has
// not been computed" (MDATK prints a KeyError as the repr of its message).

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "mdbridge/core/value.hpp"

namespace mdbridge::core {

// A numpy array: element type, shape (row-major) and values.
class NumArray {
public:
    enum class DType { Bool, Int8, Int32, Int64, Float32, Float64 };

    NumArray() = default;
    // A 1-D array, or one of the given shape (the product of `shape` must equal values.size()).
    // T: bool (stored as uint8), int8_t, int32_t, int64_t, float, double (explicitly instantiated
    // in dataset.cpp).
    template <class T>
    explicit NumArray(std::vector<T> values, std::vector<std::size_t> shape = {});

    DType dtype() const;
    std::string_view dtypeName() const;  // "bool", "int8", ... as numpy names them
    const std::vector<std::size_t>& shape() const { return m_shape; }
    std::size_t size() const;

    // The values with their own type. Throws std::invalid_argument for another type.
    template <class T>
    std::span<const T> values() const;
    // astype(float64)
    std::vector<double> toDouble() const;

    bool operator==(const NumArray&) const = default;

private:
    std::variant<std::vector<std::uint8_t>, std::vector<std::int8_t>, std::vector<std::int32_t>,
                 std::vector<std::int64_t>, std::vector<float>, std::vector<double>>
        m_data;
    std::vector<std::size_t> m_shape;
    bool m_bool = false;  // the uint8 storage holds numpy bools
};

using NamedArrays = std::vector<std::pair<std::string, NumArray>>;  // insertion order, like a dict

// One analysis' output bag (AnalysisResult).
struct AnalysisResult {
    std::string name;
    std::string key;  // cache key
    Table params;
    NamedArrays arrays;
    Table meta;  // free-form; meta["outputs"] describes the outputs

    // meta["outputs"], created as {"series": {}, "matrices": {}, "pairmaps": {}} on first use.
    Table& outputs();
    void addSeries(const std::string& out, const std::string& array, const std::string& label,
                   const std::string& unit = {}, const Table& extra = {});
    void addMatrix(const std::string& out, const std::string& array, const std::string& rows, const std::string& label,
                   const std::string& kind = "count", const Table& extra = {});
    void addPairmap(const std::string& out, const std::string& prefix, const std::string& rows, const std::string& cols,
                    const std::string& label, const Table& extra = {});

    // The array named `name`, or nullptr.
    const NumArray* array(std::string_view name) const;
    void setArray(const std::string& name, NumArray values);  // replaces, keeping the position
};

struct Series {
    std::string key;
    std::vector<double> values;
    std::string label;
    std::string unit;
    std::optional<double> target;

    // "label (unit)", or the label alone without a unit or for "count"
    std::string axisLabel() const;
};

struct Matrix {
    std::string key;
    NumArray values;  // rows x frames, its own dtype (int8 for class matrices)
    std::vector<std::string> rowLabels;
    std::string label;
    std::string kind = "count";  // count, value or class
    std::optional<std::vector<std::string>> classes;
    std::string unit;
};

// Sparse per-frame residue x residue map (CSR over frames).
struct PairMap {
    std::string key;
    std::vector<std::int64_t> offsets;  // frames + 1
    std::vector<std::int64_t> i;
    std::vector<std::int64_t> j;
    std::vector<double> values;
    std::vector<std::string> rowLabels;
    std::vector<std::string> colLabels;
    std::string label;

    std::pair<std::size_t, std::size_t> shape() const { return {rowLabels.size(), colLabels.size()}; }
    std::size_t frameCount() const { return offsets.empty() ? 0 : offsets.size() - 1; }
    // The map of frame `pos` (row-major, duplicates summed).
    std::vector<double> dense(std::size_t pos) const;
    // Fraction of frames each pair is listed in (a pair listed twice in a frame counts twice).
    std::vector<double> frequency() const;
};

// "Time (ps)", "Time (ns)", "Time (µs)", "Frame"; "Time" for other units.
std::string timeLabel(std::string_view unit);

class Dataset {
public:
    Dataset(std::vector<std::int64_t> frames, std::vector<double> timePs, std::vector<AnalysisResult> results = {},
            std::string timeUnit = "ns");

    std::size_t frameCount() const { return m_frames.size(); }
    const std::vector<std::int64_t>& frames() const { return m_frames; }
    const std::vector<double>& timePs() const { return m_timePs; }
    const std::string& timeUnit() const { return m_timeUnit; }
    const std::vector<AnalysisResult>& results() const { return m_results; }

    // Stores the result under its name (replacing one of the same name in its position).
    void add(AnalysisResult result);
    const AnalysisResult* result(std::string_view name) const;

    // The time axis in `unit` (empty: the dataset's unit; "frame": the frame numbers). Throws
    // ValueError for unknown units.
    std::vector<double> time(std::string_view unit = {}) const;
    std::string timeLabelFor(std::string_view unit = {}) const;

    // "<analysis>/<output>" keys of one kind ("series", "matrices", "pairmaps"), in insertion order.
    std::vector<std::string> keys(std::string_view kind) const;
    std::vector<std::string> seriesKeys() const { return keys("series"); }
    std::vector<std::string> matrixKeys() const { return keys("matrices"); }
    std::vector<std::string> pairmapKeys() const { return keys("pairmaps"); }
    bool has(std::string_view key) const;
    std::string kindOf(std::string_view key) const;  // "series", "matrices" or "pairmaps"

    Series series(std::string_view key) const;
    Matrix matrix(std::string_view key) const;
    PairMap pairmap(std::string_view key) const;
    // {"series": [...], "matrices": [...], "pairmaps": [...]}
    Value describe() const;

private:
    std::vector<std::int64_t> m_frames;
    std::vector<double> m_timePs;
    std::vector<AnalysisResult> m_results;
    std::string m_timeUnit;
};

// Accumulates per-frame (i, j, value) triplets into the arrays of a PairMap.
class CsrBuilder {
public:
    // One frame's entries (int32 i and j, float32 values, as MDATK stores them).
    void add(std::span<const std::int64_t> i, std::span<const std::int64_t> j, std::span<const double> v);
    // {prefix_offsets (int64), prefix_i (int32), prefix_j (int32), prefix_v (float32)}
    NamedArrays arrays(const std::string& prefix) const;

private:
    std::vector<std::int32_t> m_i, m_j;
    std::vector<float> m_v;
    std::vector<std::int64_t> m_counts;
};

}  // namespace mdbridge::core
