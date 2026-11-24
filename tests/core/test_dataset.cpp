// Analysis results and the dataset handed to panels and charts (mdatk/data.py), against what MDATK
// computes for the cases of tests/parity/make_dataset_cases.py (data/dataset, see data/README.md).

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/dataset.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

namespace {

const Array& items(const Value& value, std::string_view key) {
    return *value.find(key)->get<Array>();
}
const Array& list(const Value& value) { return *value.get<Array>(); }
const std::string& str(const Value& value) { return *value.get<std::string>(); }

double hexFloat(const Value& v) { return std::strtod(str(v).c_str(), nullptr); }

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

std::vector<double> doubles(const Array& values) {
    std::vector<double> out;
    for (const Value& v : values) {
        out.push_back(hexFloat(v));
    }
    return out;
}

std::vector<std::int64_t> ints(const Array& values) {
    std::vector<std::int64_t> out;
    for (const Value& v : values) {
        out.push_back(*v.get<std::int64_t>());
    }
    return out;
}

std::vector<std::string> strings(const Array& values) {
    std::vector<std::string> out;
    for (const Value& v : values) {
        out.push_back(str(v));
    }
    return out;
}

void expectDoubles(const std::vector<double>& got, const std::vector<double>& expected,
                   const std::string& context) {
    ASSERT_EQ(got.size(), expected.size()) << context;
    for (std::size_t k = 0; k < got.size(); ++k) {
        EXPECT_TRUE(sameDouble(got[k], expected[k]))
            << context << " " << k << ": " << got[k] << " vs " << expected[k];
    }
}

// {"dtype", "shape", "values"}
NumArray makeArray(const Value& spec) {
    const std::string& dtype = text(spec, "dtype");
    std::vector<std::size_t> shape;
    for (const Value& n : items(spec, "shape")) {
        shape.push_back(static_cast<std::size_t>(*n.get<std::int64_t>()));
    }
    const Array& values = items(spec, "values");
    if (dtype == "bool") {
        std::vector<bool> v;
        for (const Value& x : values) {
            v.push_back(*x.get<bool>());
        }
        return NumArray(std::move(v), shape);
    }
    if (dtype == "int8" || dtype == "int32" || dtype == "int64") {
        const std::vector<std::int64_t> v = ints(values);
        auto narrowed = [&v]<class T>(T) {
            std::vector<T> out;
            for (const std::int64_t x : v) {
                out.push_back(static_cast<T>(x));
            }
            return out;
        };
        if (dtype == "int8") {
            return NumArray(narrowed(std::int8_t{}), shape);
        }
        if (dtype == "int32") {
            return NumArray(narrowed(std::int32_t{}), shape);
        }
        return NumArray(v, shape);
    }
    const std::vector<double> v = doubles(values);
    if (dtype == "float32") {
        std::vector<float> f;
        for (const double x : v) {
            f.push_back(static_cast<float>(x));
        }
        return NumArray(std::move(f), shape);
    }
    return NumArray(v, shape);
}

void expectArray(const NumArray& got, const Value& expected, const std::string& context) {
    EXPECT_EQ(got.dtypeName(), text(expected, "dtype")) << context;
    std::vector<std::size_t> shape;
    for (const Value& n : items(expected, "shape")) {
        shape.push_back(static_cast<std::size_t>(*n.get<std::int64_t>()));
    }
    EXPECT_EQ(got.shape(), shape) << context;
    const Array& values = items(expected, "values");
    if (got.dtype() == NumArray::DType::Float32 || got.dtype() == NumArray::DType::Float64) {
        expectDoubles(got.toDouble(), doubles(values), context);
        return;
    }
    std::vector<std::int64_t> want;
    for (const Value& v : values) {
        const bool* b = v.get<bool>();
        want.push_back(b ? (*b ? 1 : 0) : *v.get<std::int64_t>());
    }
    std::vector<std::int64_t> have;
    for (const double v : got.toDouble()) {
        have.push_back(static_cast<std::int64_t>(v));
    }
    EXPECT_EQ(have, want) << context;
}

AnalysisResult buildResult(const Value& spec) {
    AnalysisResult result;
    result.name = text(spec, "name");
    result.key = text(spec, "key");
    if (const Table* meta = spec.find("meta")->get<Table>()) {
        result.meta = *meta;
    }
    for (const Value& entry : items(spec, "arrays")) {
        result.arrays.emplace_back(str(list(entry)[0]), makeArray(list(entry)[1]));
    }
    for (const Value& op : items(spec, "ops")) {
        const std::string& kind = str(list(op)[0]);
        const std::vector<std::string> a = strings(list(list(op)[1]));
        const Table& extra = *list(op)[2].get<Table>();
        if (kind == "add_series") {
            result.addSeries(a[0], a[1], a[2], a.size() > 3 ? a[3] : "", extra);
        } else if (kind == "add_matrix") {
            result.addMatrix(a[0], a[1], a[2], a[3], a.size() > 4 ? a[4] : "count", extra);
        } else {
            result.addPairmap(a[0], a[1], a[2], a[3], a[4], extra);
        }
    }
    return result;
}

std::string argument(const Value& arg) { return arg.isNull() ? std::string() : str(arg); }

void checkQuery(const Dataset& dataset, const std::string& op, const Value& arg,
                const Value& expected, const std::string& context) {
    const std::string a = argument(arg);
    if (const Value* error = expected.find("error")) {
        try {
            if (op == "time") {
                dataset.time(a);
            } else if (op == "kind_of") {
                dataset.kindOf(a);
            } else if (op == "series") {
                dataset.series(a);
            } else if (op == "matrix") {
                dataset.matrix(a);
            } else if (op == "pairmap") {
                dataset.pairmap(a);
            }
            ADD_FAILURE() << context << ": no error";
        } catch (const std::out_of_range& e) {
            EXPECT_EQ(text(expected, "type"), "KeyError") << context;
            EXPECT_EQ(e.what(), str(*error)) << context;
        } catch (const ValueError& e) {
            EXPECT_EQ(text(expected, "type"), "ValueError") << context;
            EXPECT_EQ(e.what(), str(*error)) << context;
        }
        return;
    }
    const Value& value = *expected.find("value");
    if (op == "time") {
        expectDoubles(dataset.time(a), doubles(list(value)), context);
    } else if (op == "time_label") {
        EXPECT_EQ(dataset.timeLabelFor(a), str(value)) << context;
    } else if (op == "keys") {
        EXPECT_EQ(dataset.keys(a), strings(list(value))) << context;
    } else if (op == "describe") {
        EXPECT_EQ(dataset.describe(), value) << context;
    } else if (op == "has") {
        EXPECT_EQ(dataset.has(a), *value.get<bool>()) << context;
    } else if (op == "kind_of") {
        EXPECT_EQ(dataset.kindOf(a), str(value)) << context;
    } else if (op == "series") {
        const Series s = dataset.series(a);
        EXPECT_EQ(s.key, text(value, "key")) << context;
        expectDoubles(s.values, doubles(items(*value.find("values"), "values")), context);
        EXPECT_EQ(s.label, text(value, "label")) << context;
        EXPECT_EQ(s.unit, text(value, "unit")) << context;
        EXPECT_EQ(s.axisLabel(), text(value, "axis_label")) << context;
        const Value& target = *value.find("target");
        ASSERT_EQ(s.target.has_value(), !target.isNull()) << context;
        if (s.target) {
            EXPECT_EQ(*s.target, hexFloat(target)) << context;
        }
    } else if (op == "matrix") {
        const Matrix m = dataset.matrix(a);
        EXPECT_EQ(m.key, text(value, "key")) << context;
        expectArray(m.values, *value.find("values"), context);
        EXPECT_EQ(m.rowLabels, strings(items(value, "row_labels"))) << context;
        EXPECT_EQ(m.label, text(value, "label")) << context;
        EXPECT_EQ(m.kind, text(value, "kind")) << context;
        EXPECT_EQ(m.unit, text(value, "unit")) << context;
        const Value& classes = *value.find("classes");
        ASSERT_EQ(m.classes.has_value(), !classes.isNull()) << context;
        if (m.classes) {
            EXPECT_EQ(*m.classes, strings(list(classes))) << context;
        }
    } else if (op == "pairmap") {
        const PairMap p = dataset.pairmap(a);
        EXPECT_EQ(p.key, text(value, "key")) << context;
        EXPECT_EQ(p.offsets, ints(items(value, "offsets"))) << context;
        EXPECT_EQ(p.i, ints(items(value, "i"))) << context;
        EXPECT_EQ(p.j, ints(items(value, "j"))) << context;
        expectDoubles(p.values, doubles(items(value, "values")), context);
        EXPECT_EQ(p.rowLabels, strings(items(value, "row_labels"))) << context;
        EXPECT_EQ(p.colLabels, strings(items(value, "col_labels"))) << context;
        EXPECT_EQ(p.label, text(value, "label")) << context;
        const Array& shape = items(value, "shape");
        EXPECT_EQ(p.shape().first, static_cast<std::size_t>(*shape[0].get<std::int64_t>()))
            << context;
        EXPECT_EQ(p.shape().second, static_cast<std::size_t>(*shape[1].get<std::int64_t>()))
            << context;
        const Array& dense = items(value, "dense");
        ASSERT_EQ(dense.size(), p.frameCount()) << context;
        for (std::size_t f = 0; f < dense.size(); ++f) {
            if (const Value* error = dense[f].find("error")) {
                EXPECT_THROW(p.dense(f), std::out_of_range) << context << " " << str(*error);
            } else {
                expectDoubles(p.dense(f), doubles(items(dense[f], "value")), context + " dense");
            }
        }
        const Value& frequency = *value.find("frequency");
        if (frequency.find("error")) {
            EXPECT_THROW(p.frequency(), std::out_of_range) << context;
        } else {
            expectDoubles(p.frequency(), doubles(items(frequency, "value")),
                          context + " frequency");
        }
    } else {
        FAIL() << op;
    }
}

}  // namespace

TEST(Dataset, MdatkCases) {
    const Value data = readJson(kData / "dataset" / "dataset_cases.json");
    int queries = 0;
    for (const Value& c : items(data, "cases")) {
        const std::string name = text(c, "name");
        std::vector<double> timePs = doubles(items(c, "time_ps"));
        Dataset dataset(ints(items(c, "frames")), timePs, {}, text(c, "time_unit"));
        for (const Value& spec : items(c, "results")) {
            AnalysisResult result = buildResult(spec);
            EXPECT_EQ(pythonRepr(Value(result.meta)), text(spec, "outputs_repr"))
                << name << " " << result.name;
            dataset.add(std::move(result));
        }
        for (const Value& q : items(c, "queries")) {
            const Array& query = list(q);
            const std::string op = str(query[0]);
            checkQuery(dataset, op, query[1], query[2], name + " " + op + " " + argument(query[1]));
            ++queries;
        }
    }
    EXPECT_GT(queries, 300);
}

TEST(Dataset, CsrBuilder) {
    const Value data = readJson(kData / "dataset" / "dataset_cases.json");
    for (const Value& c : items(data, "csr")) {
        CsrBuilder builder;
        for (const Value& frame : items(c, "frames")) {
            const Array& f = list(frame);
            const std::vector<std::int64_t> i = ints(list(f[0]));
            const std::vector<std::int64_t> j = ints(list(f[1]));
            const std::vector<double> v = doubles(list(f[2]));
            builder.add(i, j, v);
        }
        const NamedArrays arrays = builder.arrays("cm");
        const Array& expected = items(c, "arrays");
        ASSERT_EQ(arrays.size(), expected.size());
        for (std::size_t k = 0; k < arrays.size(); ++k) {
            EXPECT_EQ(arrays[k].first, str(list(expected[k])[0]));
            expectArray(arrays[k].second, list(expected[k])[1], arrays[k].first);
        }
    }
}

TEST(Dataset, AddErrors) {
    const Value data = readJson(kData / "dataset" / "dataset_cases.json");
    for (const Value& c : items(data, "add_errors")) {
        const Array& e = list(c);
        const std::string& kind = str(e[0]);
        const std::vector<std::string> a = strings(list(e[1]));
        const Table& extra = *e[2].get<Table>();
        AnalysisResult result;
        result.name = "r";
        const auto call = [&] {
            if (kind == "add_series") {
                result.addSeries(a[0], a[1], a[2], "", extra);
            } else if (kind == "add_matrix") {
                result.addMatrix(a[0], a[1], a[2], a[3], "count", extra);
            } else {
                result.addPairmap(a[0], a[1], a[2], a[3], a[4], extra);
            }
        };
        // a keyword parameter given once (unit="nm") is a C++ argument, not an extra
        if (e[3].isNull()) {
            continue;
        }
        try {
            call();
            ADD_FAILURE() << kind;
        } catch (const std::invalid_argument& error) {
            EXPECT_EQ(error.what(), str(e[3]));
        }
    }
}

TEST(Dataset, NumArray) {
    const NumArray a(std::vector<double>{1, 2, 3, 4, 5, 6}, {2, 3});
    EXPECT_EQ(a.dtype(), NumArray::DType::Float64);
    EXPECT_EQ(a.dtypeName(), "float64");
    EXPECT_EQ(a.shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_EQ(a.size(), 6u);
    EXPECT_EQ(a.values<double>()[4], 5.0);
    EXPECT_THROW(a.values<float>(), std::invalid_argument);
    EXPECT_THROW(NumArray(std::vector<double>{1, 2}, {3}), std::invalid_argument);

    const NumArray flags(std::vector<bool>{true, false});
    EXPECT_EQ(flags.dtypeName(), "bool");
    EXPECT_EQ(flags.values<std::uint8_t>()[0], 1);
    EXPECT_EQ(flags.toDouble(), (std::vector<double>{1.0, 0.0}));
    EXPECT_EQ(NumArray(std::vector<std::int8_t>{-1}).dtypeName(), "int8");
    EXPECT_EQ(NumArray(std::vector<float>{0.5f}).toDouble(), std::vector<double>{0.5});
    EXPECT_EQ(NumArray().size(), 0u);
    EXPECT_EQ(NumArray(std::vector<double>{1}).stringWidth(), 0u);
}

// numpy unicode arrays: np.array(values, dtype=str) and assignment to a '<U8' array (trace resname)
TEST(Dataset, NumArrayStrings) {
    const NumArray names(std::vector<std::string>{"G", "Mg2+", "Åx"});
    EXPECT_EQ(names.dtype(), NumArray::DType::String);
    EXPECT_EQ(names.dtypeName(), "<U4");
    EXPECT_EQ(names.stringWidth(), 4u);
    EXPECT_EQ(names.values<std::string>()[2], "Åx");
    EXPECT_THROW(names.toDouble(), std::invalid_argument);
    EXPECT_THROW(names.values<double>(), std::invalid_argument);

    const NumArray truncated = NumArray::strings({"ABCDEFGHIJ", "ééééééééé", "x中"}, 8);
    EXPECT_EQ(truncated.dtypeName(), "<U8");
    EXPECT_EQ(truncated.values<std::string>()[0], "ABCDEFGH");
    EXPECT_EQ(truncated.values<std::string>()[1], "éééééééé");
    EXPECT_EQ(truncated.values<std::string>()[2], "x中");

    const NumArray empty(std::vector<std::string>{});
    EXPECT_EQ(empty.dtypeName(), "<U1");  // np.zeros(0, dtype=str)
    EXPECT_EQ(empty.size(), 0u);
    EXPECT_EQ(NumArray(std::vector<std::string>{"", ""}).dtypeName(), "<U1");
    // numpy drops trailing NULs, the width counts them
    const NumArray nul(std::vector<std::string>{std::string("ab\0\0", 4), std::string("a\0b", 3)});
    EXPECT_EQ(nul.dtypeName(), "<U4");
    EXPECT_EQ(nul.values<std::string>()[0], "ab");
    EXPECT_EQ(nul.values<std::string>()[1], std::string("a\0b", 3));

    const NumArray grid = NumArray::strings({"a", "b", "c", "d", "e", "f"}, 2, {2, 3});
    EXPECT_EQ(grid.shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_THROW(NumArray::strings({"a"}, 1, {2}), std::invalid_argument);
    EXPECT_EQ(grid, NumArray::strings({"a", "b", "c", "d", "e", "f"}, 2, {2, 3}));
    EXPECT_NE(grid, NumArray::strings({"a", "b", "c", "d", "e", "f"}, 3, {2, 3}));
}

TEST(Dataset, ResultArraysAndOutputs) {
    AnalysisResult result;
    result.setArray("a", NumArray(std::vector<double>{1}));
    result.setArray("b", NumArray(std::vector<double>{2}));
    result.setArray("a", NumArray(std::vector<double>{3}));
    ASSERT_EQ(result.arrays.size(), 2u);
    EXPECT_EQ(result.arrays[0].first, "a");
    EXPECT_EQ(result.array("a")->values<double>()[0], 3.0);
    EXPECT_EQ(result.array("c"), nullptr);
    EXPECT_EQ(pythonRepr(Value(result.outputs())),
              "{'series': {}, 'matrices': {}, 'pairmaps': {}}");
    result.meta = {{"outputs", 3}};
    EXPECT_THROW(result.outputs(), std::invalid_argument);
}

TEST(Dataset, TimeLabels) {
    EXPECT_EQ(timeLabel("ps"), "Time (ps)");
    EXPECT_EQ(timeLabel("ns"), "Time (ns)");
    EXPECT_EQ(timeLabel("us"), "Time (µs)");
    EXPECT_EQ(timeLabel("frame"), "Frame");
    EXPECT_EQ(timeLabel("fs"), "Time");
    Series s;
    s.label = "RMSD";
    EXPECT_EQ(s.axisLabel(), "RMSD");
    s.unit = "count";
    EXPECT_EQ(s.axisLabel(), "RMSD");
    s.unit = "Å";
    EXPECT_EQ(s.axisLabel(), "RMSD (Å)");
}
