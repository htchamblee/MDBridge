// External time series (mdatk/io/external.py) and numpy.interp, against what MDATK and numpy give
// for the cases of tests/parity/make_external_cases.py (data/io/external, see data/README.md).

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/python.hpp"
#include "mdbridge/io/external.hpp"
#include "test_support.hpp"

using namespace mdbridge;
using namespace mdbridge::test;

namespace {

const fs::path kExternal = kData / "io" / "external";

const core::Array& items(const core::Value& value, std::string_view key) {
    return *value.find(key)->get<core::Array>();
}
const core::Array& list(const core::Value& value) { return *value.get<core::Array>(); }
const std::string& str(const core::Value& value) { return *value.get<std::string>(); }

double number(const core::Value& v) { return std::strtod(str(v).c_str(), nullptr); }

std::vector<double> numbers(const core::Array& values) {
    std::vector<double> out;
    for (const core::Value& v : values) {
        out.push_back(number(v));
    }
    return out;
}

std::vector<std::string> strings(const core::Array& values) {
    std::vector<std::string> out;
    for (const core::Value& v : values) {
        out.push_back(str(v));
    }
    return out;
}

bool sameDouble(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || (a == b && std::signbit(a) == std::signbit(b));
}

void expectDoubles(const std::vector<double>& got, const std::vector<double>& expected,
                   const std::string& context) {
    ASSERT_EQ(got.size(), expected.size()) << context;
    for (std::size_t k = 0; k < got.size(); ++k) {
        EXPECT_TRUE(sameDouble(got[k], expected[k]))
            << context << " [" << k << "] " << got[k] << " vs " << expected[k];
    }
}

std::string replaceAll(std::string text, std::string_view from, const std::string& to) {
    for (std::size_t at = text.find(from); at != std::string::npos;
         at = text.find(from, at + to.size())) {
        text.replace(at, from.size(), to);
    }
    return text;
}

// "{dir}/name" in MDATK's messages: the file in `dir`, with this platform's separator.
std::string withDir(const std::string& message, const TempDir& dir) {
    return replaceAll(message, "{dir}/", core::pathToUtf8(dir.path() / ""));
}

}  // namespace

TEST(External, ReadTable) {
    const core::Value cases = readJson(kExternal / "table_cases.json");
    TempDir dir;
    for (const core::Value& c : list(cases)) {
        const std::string name = text(c, "name");
        const fs::path path = dir / name;
        writeFile(path, fromHex(text(c, "bytes")));
        if (const core::Value* error = c.find("error")) {
            try {
                io::readTable(path, text(c, "delimiter"));
                ADD_FAILURE() << name << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), withDir(str(*error), dir)) << name;
            }
            continue;
        }
        const io::DataTable table = io::readTable(path, text(c, "delimiter"));
        EXPECT_EQ(table.header, strings(items(c, "header"))) << name;
        const core::Array& rows = items(c, "rows");
        ASSERT_EQ(table.rows.size(), rows.size()) << name;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            expectDoubles(table.rows[r], numbers(list(rows[r])),
                          name + " row " + std::to_string(r));
        }
    }
    EXPECT_THROW(io::readTable(dir / "missing.dat"), fs::filesystem_error);
}

TEST(External, ColumnIndexAndSanitize) {
    const core::Value cases = readJson(kExternal / "column_cases.json");
    for (const core::Value& c : items(cases, "column_index")) {
        const std::vector<std::string> header = strings(items(c, "header"));
        const std::string spec = text(c, "spec");
        const int fallback = static_cast<int>(*c.find("default")->get<std::int64_t>());
        if (const core::Value* error = c.find("error")) {
            try {
                io::columnIndex(header, spec, fallback);
                ADD_FAILURE() << spec << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), str(*error));
            }
        } else {
            EXPECT_EQ(io::columnIndex(header, spec, fallback),
                      *c.find("index")->get<std::int64_t>())
                << spec;
        }
    }
    for (const core::Value& c : items(cases, "sanitize")) {
        EXPECT_EQ(io::sanitizeKey(str(list(c)[0])), str(list(c)[1])) << str(list(c)[0]);
    }
}

TEST(External, NumpyInterp) {
    const core::Value cases = readJson(kExternal / "interp_cases.json");
    int n = 0;
    for (const core::Value& c : items(cases, "cases")) {
        const core::Array& a = list(c);
        const std::vector<double> got =
            io::numpyInterp(numbers(list(a[0])), numbers(list(a[1])), numbers(list(a[2])),
                            number(a[3]), number(a[4]));
        expectDoubles(got, numbers(list(a[5])), "case " + std::to_string(n++));
    }
    EXPECT_GT(n, 2000);
    for (const core::Value& c : items(cases, "errors")) {
        const core::Array& a = list(c);
        const auto call = [&] {
            std::vector<double> xs, ys;
            for (const core::Value& v : list(a[1])) {
                xs.push_back(*v.get<double>());
            }
            for (const core::Value& v : list(a[2])) {
                ys.push_back(*v.get<double>());
            }
            std::vector<double> t;
            for (const core::Value& v : list(a[0])) {
                t.push_back(*v.get<double>());
            }
            return io::numpyInterp(t, xs, ys, NAN, NAN);
        };
        if (a[3].isNull()) {
            EXPECT_TRUE(call().empty());
            continue;
        }
        try {
            call();
            ADD_FAILURE() << str(a[3]);
        } catch (const core::ValueError& e) {
            EXPECT_EQ(e.what(), str(a[3]));
        }
    }
}

TEST(External, LoadExternal) {
    const core::Value data = readJson(kExternal / "load_cases.json");
    TempDir dir;
    for (const auto& [name, content] : *data.find("files")->get<core::Table>()) {
        writeFile(dir / name, str(content));
    }
    for (const core::Value& c : items(data, "cases")) {
        const std::string name = text(c, "name");
        std::vector<core::ExternalSeries> specs;
        for (const core::Value& s : items(c, "specs")) {
            core::ExternalSeries spec;
            const std::string& path = text(s, "path");
            spec.path = path.empty() ? "" : core::pathToUtf8(dir.path() / path);
            spec.label = text(s, "label");
            spec.x_column = text(s, "x_column");
            if (const core::Value* unit = s.find("x_unit")) {
                spec.x_unit = str(*unit);
            }
            if (const core::Value* columns = s.find("columns")) {
                spec.columns = strings(list(*columns));
            }
            spec.delimiter = text(s, "delimiter");
            specs.push_back(spec);
        }
        std::vector<std::int64_t> frames;
        for (const core::Value& f : items(c, "frames")) {
            frames.push_back(*f.get<std::int64_t>());
        }
        const core::Dataset dataset(frames, numbers(items(c, "time_ps")));
        if (const core::Value* error = c.find("error")) {
            if (text(c, "error_type") == "FileNotFoundError") {
                EXPECT_THROW(io::loadExternal(specs, dataset), fs::filesystem_error) << name;
                continue;
            }
            try {
                io::loadExternal(specs, dataset);
                ADD_FAILURE() << name << ": no error";
            } catch (const core::ValueError& e) {
                EXPECT_EQ(e.what(), withDir(str(*error), dir)) << name;
            }
            continue;
        }
        const core::AnalysisResult result = io::loadExternal(specs, dataset);
        EXPECT_EQ(result.name, "external");
        const core::Array& arrays = items(c, "arrays");
        ASSERT_EQ(result.arrays.size(), arrays.size()) << name;
        for (std::size_t k = 0; k < arrays.size(); ++k) {
            EXPECT_EQ(result.arrays[k].first, str(list(arrays[k])[0])) << name;
            EXPECT_EQ(result.arrays[k].second.dtypeName(), "float64") << name;
            expectDoubles(result.arrays[k].second.toDouble(), numbers(list(list(arrays[k])[1])),
                          name);
        }
        EXPECT_EQ(core::pythonRepr(core::Value(result.meta)), text(c, "meta_repr")) << name;
    }
}

// MDATK's tests/test_render_cache.py::test_external_table_and_interpolation
TEST(External, MdatkTableAndInterpolation) {
    TempDir dir;
    writeFile(dir / "rmsd.dat", "#Frame  rmsd  other\n1 0.0 5\n2 1.0 6\n3 2.0 7\n");
    const io::DataTable table = io::readTable(dir / "rmsd.dat");
    EXPECT_EQ(table.header, (std::vector<std::string>{"Frame", "rmsd", "other"}));
    ASSERT_EQ(table.rows.size(), 3u);
    EXPECT_EQ(table.rows[0].size(), 3u);

    core::Dataset dataset({0, 1, 2}, {0.0, 10.0, 20.0});
    core::ExternalSeries spec;
    spec.path = (dir.path() / "rmsd.dat").string();
    spec.label = "cpp";
    spec.x_column = "Frame";
    spec.x_unit = "frame1";
    spec.columns = {"rmsd"};
    dataset.add(io::loadExternal({spec}, dataset));
    EXPECT_EQ(dataset.series("external/cpp.rmsd").values, (std::vector<double>{0, 1, 2}));

    writeFile(dir / "e.csv", "time_ns,energy\n0,1\n0.02,3\n");
    core::ExternalSeries csv;
    csv.path = (dir.path() / "e.csv").string();
    csv.x_unit = "ns";
    const core::AnalysisResult result = io::loadExternal({csv}, dataset);
    const std::vector<double> values = result.arrays.front().second.toDouble();
    ASSERT_EQ(values.size(), 3u);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_NEAR(values[k], static_cast<double>(k + 1), 1e-7);
    }
}
