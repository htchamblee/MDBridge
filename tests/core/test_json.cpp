#include <limits>

#include <gtest/gtest.h>

#include "mdbridge/core/json.hpp"

using namespace mdbridge::core;

TEST(PythonFloatRepr, MatchesPythonRepr) {
    // expected values: repr(x) in CPython
    const std::pair<double, const char*> cases[] = {
        {0.0, "0.0"},
        {-0.0, "-0.0"},
        {1.0, "1.0"},
        {2.5, "2.5"},
        {-3.25, "-3.25"},
        {0.1, "0.1"},
        {1.0 / 3.0, "0.3333333333333333"},
        {100000.0, "100000.0"},
        {123456.789, "123456.789"},
        {1e15, "1000000000000000.0"},
        {1e16, "1e+16"},
        {1.5e16, "1.5e+16"},
        {0.0001, "0.0001"},
        {0.00012, "0.00012"},
        {1e-05, "1e-05"},
        {1e-07, "1e-07"},
        {1.25e-10, "1.25e-10"},
        {1e100, "1e+100"},
        {5e-324, "5e-324"},
        {std::numeric_limits<double>::max(), "1.7976931348623157e+308"},
        {std::numeric_limits<double>::quiet_NaN(), "nan"},
        {std::numeric_limits<double>::infinity(), "inf"},
        {-std::numeric_limits<double>::infinity(), "-inf"},
    };
    for (const auto& [value, expected] : cases) {
        EXPECT_EQ(pythonFloatRepr(value), expected) << expected;
    }
}

TEST(PythonStringRepr, ChoosesQuotesAndEscapes) {
    EXPECT_EQ(pythonStringRepr("something-else"), "'something-else'");
    EXPECT_EQ(pythonStringRepr("it's"), "\"it's\"");
    EXPECT_EQ(pythonStringRepr("it's \"x\""), "'it\\'s \"x\"'");
    EXPECT_EQ(pythonStringRepr("a\\b\n\t\x01"), "'a\\\\b\\n\\t\\x01'");
    EXPECT_EQ(pythonStringRepr("Mg²⁺"), "'Mg²⁺'");
}

TEST(ToJson, FormatsLikePythonJsonDumps) {
    const Value value = Table{
        {"s", "quote \" backslash \\ tab\t newline\n control \x01 Å"},
        {"i", -3},
        {"f", 2.5},
        {"b", true},
        {"none", Value()},
        {"empty_list", Array{}},
        {"empty_table", Table{}},
        {"list", Array{1, 2.0, "x"}},
        {"nested", Table{{"nan", std::numeric_limits<double>::quiet_NaN()},
                         {"inf", -std::numeric_limits<double>::infinity()}}},
    };
    EXPECT_EQ(toJson(value),
              "{\n"
              "  \"s\": \"quote \\\" backslash \\\\ tab\\t newline\\n control \\u0001 Å\",\n"
              "  \"i\": -3,\n"
              "  \"f\": 2.5,\n"
              "  \"b\": true,\n"
              "  \"none\": null,\n"
              "  \"empty_list\": [],\n"
              "  \"empty_table\": {},\n"
              "  \"list\": [\n"
              "    1,\n"
              "    2.0,\n"
              "    \"x\"\n"
              "  ],\n"
              "  \"nested\": {\n"
              "    \"nan\": NaN,\n"
              "    \"inf\": -Infinity\n"
              "  }\n"
              "}");
    EXPECT_EQ(toJson(Table{{"a", Array{1, Table{{"b", false}}}}}, -1), "{\"a\": [1, {\"b\": false}]}");
}

TEST(DateTime, IsoformatLikePython) {
    DateTime dt{DateTime::Kind::LocalDateTime, 2026, 10, 15, 14, 32, 5, 0, 0};
    EXPECT_EQ(dt.isoformat(), "2026-10-15T14:32:05");
    dt.microsecond = 250000;
    EXPECT_EQ(dt.isoformat(), "2026-10-15T14:32:05.250000");
    dt.kind = DateTime::Kind::OffsetDateTime;
    dt.offsetMinutes = -330;
    EXPECT_EQ(dt.isoformat(), "2026-10-15T14:32:05.250000-05:30");
    dt.offsetMinutes = 0;
    EXPECT_EQ(dt.isoformat(), "2026-10-15T14:32:05.250000+00:00");
    EXPECT_EQ((DateTime{DateTime::Kind::LocalDate, 2026, 1, 2}).isoformat(), "2026-01-02");
    EXPECT_EQ((DateTime{DateTime::Kind::LocalTime, 0, 0, 0, 7, 3, 0}).isoformat(), "07:03:00");
}
