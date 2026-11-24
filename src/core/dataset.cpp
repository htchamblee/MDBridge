#include "mdbridge/core/dataset.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/core/selections.hpp"
#include "mdbridge/core/text.hpp"

namespace mdbridge::core {

namespace {

// d[key] = value of a dict: replaced in its position, or appended
void assign(Table& table, const std::string& key, Value value) {
    for (auto& [k, v] : table) {
        if (k == key) {
            v = std::move(value);
            return;
        }
    }
    table.emplace_back(key, std::move(value));
}

const Value* lookup(const Table& table, std::string_view key) {
    for (const auto& [k, v] : table) {
        if (k == key) {
            return &v;
        }
    }
    return nullptr;
}

// d.setdefault(key, {}) for a dict value
Table& tableEntry(Table& table, const std::string& key) {
    for (auto& [k, v] : table) {
        if (k == key) {
            if (Table* t = v.get<Table>()) {
                return *t;
            }
            throw std::invalid_argument("'" + key + "' is not a table");
        }
    }
    table.emplace_back(key, Table{});
    return *table.back().second.get<Table>();
}

// dict(array=..., label=..., **extra) of a call add_x(self, out, <fields>, **extra): an extra
// keyword naming a parameter is the TypeError of the call
Table descriptor(std::string_view method, Table fields, const Table& extra) {
    for (const auto& [k, v] : extra) {
        if (k == "self" || k == "out" || lookup(fields, k)) {
            throw std::invalid_argument("AnalysisResult." + std::string(method) +
                                        "() got multiple values for argument '" + k + "'");
        }
        fields.emplace_back(k, v);
    }
    return fields;
}

std::string textOf(const Value& value) {
    if (const std::string* s = value.get<std::string>()) {
        return *s;
    }
    return pythonStr(value);
}

// d.get(key, fallback) as text
std::string textOr(const Table& d, std::string_view key, std::string_view fallback) {
    const Value* v = lookup(d, key);
    return v ? textOf(*v) : std::string(fallback);
}

// d[key] (KeyError)
const Value& required(const Table& d, std::string_view key) {
    const Value* v = lookup(d, key);
    if (!v) {
        throw std::out_of_range(std::string(key));
    }
    return *v;
}

// list(value) of a meta entry: the items of a list, the characters of a string
std::vector<std::string> labelList(const Value* value) {
    std::vector<std::string> out;
    if (!value) {
        return out;
    }
    if (const Array* a = value->get<Array>()) {
        for (const Value& item : *a) {
            out.push_back(textOf(item));
        }
    } else if (const std::string* s = value->get<std::string>()) {
        for (const char32_t c : decodeUtf8(*s)) {
            out.push_back(encodeUtf8(std::u32string(1, c)));
        }
    } else if (const Table* t = value->get<Table>()) {
        for (const auto& [k, v] : *t) {
            out.push_back(k);
        }
    }
    return out;
}

const Table* outputsOf(const AnalysisResult& result) {
    const Value* outputs = lookup(result.meta, "outputs");
    return outputs ? outputs->get<Table>() : nullptr;
}

const NumArray& arrayOf(const AnalysisResult& result, const std::string& name) {
    const NumArray* array = result.array(name);
    if (!array) {
        throw std::out_of_range(name);
    }
    return *array;
}

std::vector<std::int64_t> toInt64(const NumArray& array) {
    std::vector<std::int64_t> out;
    out.reserve(array.size());
    switch (array.dtype()) {
    case NumArray::DType::Bool:
        for (const std::uint8_t v : array.values<std::uint8_t>()) {
            out.push_back(v);
        }
        break;
    case NumArray::DType::Int8:
        for (const std::int8_t v : array.values<std::int8_t>()) {
            out.push_back(v);
        }
        break;
    case NumArray::DType::Int32:
        for (const std::int32_t v : array.values<std::int32_t>()) {
            out.push_back(v);
        }
        break;
    case NumArray::DType::Int64: {
        const auto v = array.values<std::int64_t>();
        out.assign(v.begin(), v.end());
        break;
    }
    default:
        for (const double v : array.toDouble()) {
            out.push_back(static_cast<std::int64_t>(v));
        }
    }
    return out;
}

std::size_t numpyIndex(std::int64_t index, std::size_t size) {
    const auto n = static_cast<std::int64_t>(size);
    if (index < -n || index >= n) {
        throw std::out_of_range("index " + std::to_string(index) +
                                " is out of bounds for axis 0 with size " + std::to_string(size));
    }
    return static_cast<std::size_t>(index < 0 ? index + n : index);
}

}  // namespace

// ---------------------------------------------------------------- NumArray

template <class T>
NumArray::NumArray(std::vector<T> values, std::vector<std::size_t> shape)
    : m_shape(std::move(shape)) {
    const std::size_t count = values.size();
    if constexpr (std::is_same_v<T, bool>) {
        m_data = std::vector<std::uint8_t>(values.begin(), values.end());
        m_bool = true;
    } else {
        m_data = std::move(values);
    }
    if (m_shape.empty()) {
        m_shape = {count};
    }
    std::size_t product = 1;
    for (const std::size_t n : m_shape) {
        product *= n;
    }
    if (product != count) {
        throw std::invalid_argument("NumArray: the shape does not match the number of values");
    }
}

template NumArray::NumArray(std::vector<bool>, std::vector<std::size_t>);
template NumArray::NumArray(std::vector<std::int8_t>, std::vector<std::size_t>);
template NumArray::NumArray(std::vector<std::int32_t>, std::vector<std::size_t>);
template NumArray::NumArray(std::vector<std::int64_t>, std::vector<std::size_t>);
template NumArray::NumArray(std::vector<float>, std::vector<std::size_t>);
template NumArray::NumArray(std::vector<double>, std::vector<std::size_t>);

NumArray::DType NumArray::dtype() const {
    switch (m_data.index()) {
    case 0: return DType::Bool;  // the uint8 storage only holds bools
    case 1: return DType::Int8;
    case 2: return DType::Int32;
    case 3: return DType::Int64;
    case 4: return DType::Float32;
    default: return DType::Float64;
    }
}

std::string_view NumArray::dtypeName() const {
    switch (dtype()) {
    case DType::Bool: return "bool";
    case DType::Int8: return "int8";
    case DType::Int32: return "int32";
    case DType::Int64: return "int64";
    case DType::Float32: return "float32";
    case DType::Float64: return "float64";
    }
    return "";
}

std::size_t NumArray::size() const {
    return std::visit([](const auto& v) { return v.size(); }, m_data);
}

template <class T>
std::span<const T> NumArray::values() const {
    if (const auto* v = std::get_if<std::vector<T>>(&m_data)) {
        return *v;
    }
    throw std::invalid_argument("NumArray: the values are " + std::string(dtypeName()));
}

template std::span<const std::uint8_t> NumArray::values() const;
template std::span<const std::int8_t> NumArray::values() const;
template std::span<const std::int32_t> NumArray::values() const;
template std::span<const std::int64_t> NumArray::values() const;
template std::span<const float> NumArray::values() const;
template std::span<const double> NumArray::values() const;

std::vector<double> NumArray::toDouble() const {
    return std::visit(
        [](const auto& v) {
            std::vector<double> out;
            out.reserve(v.size());
            for (const auto x : v) {
                out.push_back(static_cast<double>(x));
            }
            return out;
        },
        m_data);
}

// ---------------------------------------------------------------- AnalysisResult

Table& AnalysisResult::outputs() {
    for (auto& [k, v] : meta) {
        if (k == "outputs") {
            if (Table* t = v.get<Table>()) {
                return *t;
            }
            throw std::invalid_argument("meta['outputs'] is not a table");
        }
    }
    meta.emplace_back("outputs",
                      Table{{"series", Table{}}, {"matrices", Table{}}, {"pairmaps", Table{}}});
    return *meta.back().second.get<Table>();
}

void AnalysisResult::addSeries(const std::string& out, const std::string& array,
                               const std::string& label, const std::string& unit,
                               const Table& extra) {
    Table d = descriptor("add_series", {{"array", array}, {"label", label}, {"unit", unit}}, extra);
    assign(tableEntry(outputs(), "series"), out, std::move(d));
}

void AnalysisResult::addMatrix(const std::string& out, const std::string& array,
                               const std::string& rows, const std::string& label,
                               const std::string& kind, const Table& extra) {
    Table d = descriptor(
        "add_matrix", {{"array", array}, {"rows", rows}, {"label", label}, {"kind", kind}}, extra);
    assign(tableEntry(outputs(), "matrices"), out, std::move(d));
}

void AnalysisResult::addPairmap(const std::string& out, const std::string& prefix,
                                const std::string& rows, const std::string& cols,
                                const std::string& label, const Table& extra) {
    Table d =
        descriptor("add_pairmap",
                   {{"prefix", prefix}, {"rows", rows}, {"cols", cols}, {"label", label}}, extra);
    assign(tableEntry(outputs(), "pairmaps"), out, std::move(d));
}

const NumArray* AnalysisResult::array(std::string_view name) const {
    for (const auto& [k, v] : arrays) {
        if (k == name) {
            return &v;
        }
    }
    return nullptr;
}

void AnalysisResult::setArray(const std::string& name, NumArray values) {
    for (auto& [k, v] : arrays) {
        if (k == name) {
            v = std::move(values);
            return;
        }
    }
    arrays.emplace_back(name, std::move(values));
}

// ---------------------------------------------------------------- Series, PairMap

std::string Series::axisLabel() const {
    return !unit.empty() && unit != "count" ? label + " (" + unit + ")" : label;
}

std::vector<double> PairMap::dense(std::size_t pos) const {
    const auto [rows, cols] = shape();
    std::vector<double> out(rows * cols, 0.0);
    if (pos + 1 >= offsets.size()) {
        throw std::out_of_range("index " + std::to_string(pos + 1) +
                                " is out of bounds for axis 0 with size " +
                                std::to_string(offsets.size()));
    }
    const auto clamp = [&](std::int64_t x, std::size_t n) {
        const auto size = static_cast<std::int64_t>(n);
        if (x < 0) {
            x = std::max<std::int64_t>(x + size, 0);
        }
        return static_cast<std::size_t>(std::min(x, size));
    };
    const std::size_t s = clamp(offsets[pos], i.size());
    const std::size_t e = clamp(offsets[pos + 1], i.size());
    for (std::size_t k = s; k < e; ++k) {
        if (k >= j.size() || k >= values.size()) {
            throw std::invalid_argument("PairMap: i, j and values differ in length");
        }
        out[numpyIndex(i[k], rows) * cols + numpyIndex(j[k], cols)] += values[k];
    }
    return out;
}

std::vector<double> PairMap::frequency() const {
    const auto [rows, cols] = shape();
    std::vector<double> out(rows * cols, 0.0);
    const auto frames = static_cast<std::int64_t>(offsets.size()) - 1;
    if (frames == 0) {
        return out;
    }
    if (i.size() != j.size()) {
        throw std::invalid_argument("PairMap: i and j differ in length");
    }
    for (std::size_t k = 0; k < i.size(); ++k) {
        out[numpyIndex(i[k], rows) * cols + numpyIndex(j[k], cols)] += 1.0;
    }
    for (double& v : out) {
        v /= static_cast<double>(frames);
    }
    return out;
}

// ---------------------------------------------------------------- Dataset

std::string timeLabel(std::string_view unit) {
    if (unit == "ps") {
        return "Time (ps)";
    }
    if (unit == "ns") {
        return "Time (ns)";
    }
    if (unit == "us") {
        return "Time (µs)";
    }
    if (unit == "frame") {
        return "Frame";
    }
    return "Time";
}

Dataset::Dataset(std::vector<std::int64_t> frames, std::vector<double> timePs,
                 std::vector<AnalysisResult> results, std::string timeUnit)
    : m_frames(std::move(frames)), m_timePs(std::move(timePs)), m_timeUnit(std::move(timeUnit)) {
    for (AnalysisResult& result : results) {
        add(std::move(result));
    }
}

void Dataset::add(AnalysisResult result) {
    for (AnalysisResult& r : m_results) {
        if (r.name == result.name) {
            r = std::move(result);
            return;
        }
    }
    m_results.push_back(std::move(result));
}

const AnalysisResult* Dataset::result(std::string_view name) const {
    for (const AnalysisResult& r : m_results) {
        if (r.name == name) {
            return &r;
        }
    }
    return nullptr;
}

std::vector<double> Dataset::time(std::string_view unit) const {
    const std::string_view u = unit.empty() ? std::string_view(m_timeUnit) : unit;
    std::vector<double> out;
    if (u == "frame") {
        for (const std::int64_t f : m_frames) {
            out.push_back(static_cast<double>(f));
        }
        return out;
    }
    const double factor = timeUnitFactor(u);
    for (const double t : m_timePs) {
        out.push_back(t / factor);
    }
    return out;
}

std::string Dataset::timeLabelFor(std::string_view unit) const {
    return timeLabel(unit.empty() ? std::string_view(m_timeUnit) : unit);
}

namespace {

struct Found {
    const AnalysisResult* result;
    const Table* desc;
};

Found lookupKey(const Dataset& dataset, std::string_view key, std::string_view kind) {
    const std::size_t slash = key.find('/');
    if (slash == std::string_view::npos) {
        throw std::out_of_range("data key " + pythonStringRepr(std::string(key)) +
                                " must look like '<analysis>/<output>'");
    }
    const std::string name(key.substr(0, slash));
    const std::string out(key.substr(slash + 1));
    const AnalysisResult* result = dataset.result(name);
    if (!result) {
        throw std::out_of_range("analysis " + pythonStringRepr(name) + " has not been computed");
    }
    const Table* outputs = outputsOf(*result);
    const Value* group = outputs ? lookup(*outputs, kind) : nullptr;
    const Value* desc = group && group->isTable() ? group->find(out) : nullptr;
    if (!desc || desc->isNull()) {
        throw std::out_of_range(name + " has no " + std::string(kind.substr(0, kind.size() - 1)) +
                                " output " + pythonStringRepr(out));
    }
    const Table* table = desc->get<Table>();
    if (!table) {
        throw std::invalid_argument("the description of " + std::string(key) + " is not a table");
    }
    return {result, table};
}

}  // namespace

std::vector<std::string> Dataset::keys(std::string_view kind) const {
    std::vector<std::string> out;
    for (const AnalysisResult& result : m_results) {
        const Table* outputs = outputsOf(result);
        const Value* group = outputs ? lookup(*outputs, kind) : nullptr;
        if (const Table* t = group ? group->get<Table>() : nullptr) {
            for (const auto& [k, v] : *t) {
                out.push_back(result.name + "/" + k);
            }
        }
    }
    return out;
}

bool Dataset::has(std::string_view key) const {
    try {
        kindOf(key);
        return true;
    } catch (const std::out_of_range&) {
        return false;
    }
}

std::string Dataset::kindOf(std::string_view key) const {
    for (const char* kind : {"series", "matrices", "pairmaps"}) {
        try {
            lookupKey(*this, key, kind);
            return kind;
        } catch (const std::out_of_range&) {
        }
    }
    throw std::out_of_range(std::string(key));
}

Series Dataset::series(std::string_view key) const {
    const auto [result, d] = lookupKey(*this, key, "series");
    Series s;
    s.key = key;
    s.values = arrayOf(*result, textOf(required(*d, "array"))).toDouble();
    s.label = textOr(*d, "label", key);
    s.unit = textOr(*d, "unit", "");
    if (const Value* target = lookup(*d, "target")) {
        if (const double* x = target->get<double>()) {
            s.target = *x;
        } else if (const std::int64_t* n = target->get<std::int64_t>()) {
            s.target = static_cast<double>(*n);
        } else if (const bool* b = target->get<bool>()) {
            s.target = *b ? 1.0 : 0.0;
        }
    }
    return s;
}

Matrix Dataset::matrix(std::string_view key) const {
    const auto [result, d] = lookupKey(*this, key, "matrices");
    Matrix m;
    m.key = key;
    m.values = arrayOf(*result, textOf(required(*d, "array")));
    m.rowLabels = labelList(lookup(result->meta, textOf(required(*d, "rows"))));
    m.label = textOr(*d, "label", key);
    m.kind = textOr(*d, "kind", "count");
    if (const Value* classes = lookup(*d, "classes"); classes && !classes->isNull()) {
        m.classes = labelList(classes);
    }
    m.unit = textOr(*d, "unit", "");
    return m;
}

PairMap Dataset::pairmap(std::string_view key) const {
    const auto [result, d] = lookupKey(*this, key, "pairmaps");
    const std::string p = textOf(required(*d, "prefix"));
    PairMap map;
    map.key = key;
    map.offsets = toInt64(arrayOf(*result, p + "_offsets"));
    map.i = toInt64(arrayOf(*result, p + "_i"));
    map.j = toInt64(arrayOf(*result, p + "_j"));
    map.values = arrayOf(*result, p + "_v").toDouble();
    map.rowLabels = labelList(lookup(result->meta, textOf(required(*d, "rows"))));
    map.colLabels = labelList(lookup(result->meta, textOf(required(*d, "cols"))));
    map.label = textOr(*d, "label", key);
    return map;
}

Value Dataset::describe() const {
    const auto list = [](const std::vector<std::string>& keys) {
        Array out;
        for (const std::string& k : keys) {
            out.emplace_back(k);
        }
        return out;
    };
    return Table{{"series", list(seriesKeys())},
                 {"matrices", list(matrixKeys())},
                 {"pairmaps", list(pairmapKeys())}};
}

// ---------------------------------------------------------------- CsrBuilder

void CsrBuilder::add(std::span<const std::int64_t> i, std::span<const std::int64_t> j,
                     std::span<const double> v) {
    for (const std::int64_t x : i) {
        m_i.push_back(static_cast<std::int32_t>(x));
    }
    for (const std::int64_t x : j) {
        m_j.push_back(static_cast<std::int32_t>(x));
    }
    for (const double x : v) {
        m_v.push_back(static_cast<float>(x));
    }
    m_counts.push_back(static_cast<std::int64_t>(i.size()));
}

NamedArrays CsrBuilder::arrays(const std::string& prefix) const {
    std::vector<std::int64_t> offsets(m_counts.size() + 1, 0);
    for (std::size_t k = 0; k < m_counts.size(); ++k) {
        offsets[k + 1] = offsets[k] + m_counts[k];
    }
    return {{prefix + "_offsets", NumArray(std::move(offsets))},
            {prefix + "_i", NumArray(m_i)},
            {prefix + "_j", NumArray(m_j)},
            {prefix + "_v", NumArray(m_v)}};
}

}  // namespace mdbridge::core
