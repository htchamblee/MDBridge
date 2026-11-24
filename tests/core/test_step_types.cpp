// Step types of the Project Planner (MDATK's mdatk/workflows schema.py, base.py, registry.py and
// the settings comments of projectfile.py). Expected values made with MDATK v0.0.4.1 by
// tools/gen_step_types.py and tests/parity/make_step_type_cases.py (see data/README.md).
//
// The registry is process-wide. gtest runs the suites in the order they first appear and the tests
// of a suite in order: StepTypeRegistry lists the built-in step types first, the comment tests and
// StepTypeRegister (last) add step types.

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/python.hpp"
#include "mdbridge/core/step_types.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using mdbridge::test::kData;
using mdbridge::test::readJson;

namespace {

Value cases(const std::string& name) { return readJson(kData / "step-types" / name); }

const Array& items(const Value& value, std::string_view key) {
    return *value.find(key)->get<Array>();
}
std::string str(const Value& value, std::string_view key) { return value.find(key)->stringOr(); }
bool flag(const Value& value, std::string_view key) {
    const Value* v = value.find(key);
    return v && v->get<bool>() && *v->get<bool>();
}
std::optional<std::string> optionalString(const Value& value) {
    if (value.isNull()) {
        return std::nullopt;
    }
    return *value.get<std::string>();
}
std::optional<double> number(const Value& value) {
    if (const double* f = value.get<double>()) {
        return *f;
    }
    if (const std::int64_t* i = value.get<std::int64_t>()) {
        return static_cast<double>(*i);
    }
    return std::nullopt;
}
std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> out;
    for (const Value& v : *value.get<Array>()) {
        out.push_back(*v.get<std::string>());
    }
    return out;
}

// The tagged values of make_step_type_cases.py
Value decode(const Value& v) {
    const std::string t = str(v, "t");
    if (t == "none") {
        return {};
    }
    if (t == "bool") {
        return *v.find("v");
    }
    if (t == "int") {
        return Value(std::int64_t{std::stoll(str(v, "v"))});
    }
    if (t == "float") {
        const std::string hex = str(v, "hex");
        if (hex == "nan") {
            return std::nan("");
        }
        if (hex == "inf" || hex == "-inf") {
            return (hex[0] == '-' ? -1 : 1) * HUGE_VAL;
        }
        return std::strtod(hex.c_str(), nullptr);
    }
    if (t == "str") {
        return str(v, "v");
    }
    if (t == "datetime") {
        return mdbridge::test::parseIso(str(v, "iso"), DateTime::Kind::LocalDateTime);
    }
    if (t == "date") {
        return mdbridge::test::parseIso(str(v, "iso"), DateTime::Kind::LocalDate);
    }
    if (t == "time") {
        return mdbridge::test::parseIso(str(v, "iso"), DateTime::Kind::LocalTime);
    }
    if (t == "array") {
        Array out;
        for (const Value& item : items(v, "v")) {
            out.push_back(decode(item));
        }
        return out;
    }
    if (t == "table") {
        Table out;
        for (const Value& pair : items(v, "v")) {
            const Array& kv = *pair.get<Array>();
            out.emplace_back(*kv[0].get<std::string>(), decode(kv[1]));
        }
        return out;
    }
    throw std::runtime_error("unknown tag " + t);
}

// Value equality that also tells NaN from NaN and 0.0 from -0.0
bool same(const Value& a, const Value& b) {
    const double* fa = a.get<double>();
    const double* fb = b.get<double>();
    if (fa && fb) {
        if (std::isnan(*fa) || std::isnan(*fb)) {
            return std::isnan(*fa) && std::isnan(*fb);
        }
        return *fa == *fb && std::signbit(*fa) == std::signbit(*fb);
    }
    if (a.type() != b.type()) {
        return false;
    }
    if (const Array* x = a.get<Array>()) {
        const Array& y = *b.get<Array>();
        if (x->size() != y.size()) {
            return false;
        }
        for (std::size_t i = 0; i < x->size(); ++i) {
            if (!same((*x)[i], y[i])) {
                return false;
            }
        }
        return true;
    }
    if (const Table* x = a.get<Table>()) {
        const Table& y = *b.get<Table>();
        if (x->size() != y.size()) {
            return false;
        }
        for (std::size_t i = 0; i < x->size(); ++i) {
            if ((*x)[i].first != y[i].first || !same((*x)[i].second, y[i].second)) {
                return false;
            }
        }
        return true;
    }
    return a == b;
}

Requirements requirements(const Value& pairs) {
    Requirements out;
    for (const Value& pair : *pairs.get<Array>()) {
        const Array& kv = *pair.get<Array>();
        out.emplace_back(*kv[0].get<std::string>(), decode(kv[1]));
    }
    return out;
}

OptionKind kindOf(const std::string& name) {
    for (OptionKind k : {OptionKind::Bool, OptionKind::Choice, OptionKind::Radio, OptionKind::Int,
                         OptionKind::Float, OptionKind::Text, OptionKind::Multiline,
                         OptionKind::File}) {
        if (optionKindName(k) == name) {
            return k;
        }
    }
    throw std::runtime_error("unknown kind " + name);
}

Option optionFrom(const Value& j) {
    Option o;
    o.key = str(j, "key");
    o.label = str(j, "label");
    o.kind = kindOf(str(j, "kind"));
    o.defaultValue = decode(*j.find("default"));
    o.tooltip = str(j, "tooltip");
    for (const Value& c : items(j, "choices")) {
        const Array& pair = *c.get<Array>();
        o.choices.push_back({*pair[0].get<std::string>(), *pair[1].get<std::string>()});
    }
    o.minimum = number(*j.find("minimum"));
    o.maximum = number(*j.find("maximum"));
    o.step = number(*j.find("step"));
    o.digits = static_cast<int>(*j.find("digits")->get<std::int64_t>());
    o.optional = flag(j, "optional");
    o.placeholder = str(j, "placeholder");
    o.unit = str(j, "unit");
    o.requirements = requirements(*j.find("requires"));
    return o;
}

Port portFrom(const Value& j) {
    return Port{str(j, "name"), str(j, "label"), strings(*j.find("kinds")), flag(j, "optional"),
                flag(j, "multiple"), str(j, "tooltip")};
}

StepType stepTypeFrom(const Value& j) {
    StepType t;
    t.name = str(j, "name");
    t.label = str(j, "label");
    t.category = str(j, "category");
    t.tag = str(j, "tag");
    t.description = str(j, "description");
    t.generates = flag(j, "generates");
    t.pathSettings = strings(*j.find("path_settings"));
    for (const Value& p : items(j, "inputs")) {
        t.inputs.push_back(portFrom(p));
    }
    for (const Value& p : items(j, "outputs")) {
        t.outputs.push_back(portFrom(p));
    }
    for (const Value& sj : items(j, "sections")) {
        Section s;
        s.key = str(sj, "key");
        s.title = str(sj, "title");
        s.tooltip = str(sj, "tooltip");
        s.repeat = flag(sj, "repeat");
        s.itemLabel = str(sj, "item_label");
        for (const Value& o : items(sj, "options")) {
            s.options.push_back(optionFrom(o));
        }
        for (const Value& tj : items(sj, "templates")) {
            s.templates.push_back(
                {str(tj, "id"), str(tj, "label"), *decode(*tj.find("values")).get<Table>()});
        }
        s.defaultItems = strings(*sj.find("default_items"));
        s.requirements = requirements(*sj.find("requires"));
        t.sections.push_back(std::move(s));
    }
    return t;
}

void expectSameOption(const Option& a, const Option& b) {
    SCOPED_TRACE("option " + b.key);
    EXPECT_EQ(a.key, b.key);
    EXPECT_EQ(a.label, b.label);
    EXPECT_EQ(a.kind, b.kind);
    EXPECT_TRUE(same(a.defaultValue, b.defaultValue)) << pythonRepr(a.defaultValue);
    EXPECT_EQ(a.tooltip, b.tooltip);
    ASSERT_EQ(a.choices.size(), b.choices.size());
    for (std::size_t i = 0; i < a.choices.size(); ++i) {
        EXPECT_EQ(a.choices[i].value, b.choices[i].value);
        EXPECT_EQ(a.choices[i].label, b.choices[i].label);
    }
    EXPECT_EQ(a.minimum, b.minimum);
    EXPECT_EQ(a.maximum, b.maximum);
    EXPECT_EQ(a.step, b.step);
    EXPECT_EQ(a.digits, b.digits);
    EXPECT_EQ(a.optional, b.optional);
    EXPECT_EQ(a.placeholder, b.placeholder);
    EXPECT_EQ(a.unit, b.unit);
    EXPECT_EQ(a.requirements, b.requirements);
}

void expectSamePort(const Port& a, const Port& b) {
    SCOPED_TRACE("port " + b.name);
    EXPECT_EQ(a.name, b.name);
    EXPECT_EQ(a.label, b.label);
    EXPECT_EQ(a.kinds, b.kinds);
    EXPECT_EQ(a.optional, b.optional);
    EXPECT_EQ(a.multiple, b.multiple);
    EXPECT_EQ(a.tooltip, b.tooltip);
}

void expectSameStepType(const StepType& a, const StepType& b) {
    SCOPED_TRACE("step type " + b.name);
    EXPECT_EQ(a.name, b.name);
    EXPECT_EQ(a.label, b.label);
    EXPECT_EQ(a.category, b.category);
    EXPECT_EQ(a.tag, b.tag);
    EXPECT_EQ(a.description, b.description);
    EXPECT_EQ(a.generates, b.generates);
    EXPECT_EQ(a.pathSettings, b.pathSettings);
    ASSERT_EQ(a.inputs.size(), b.inputs.size());
    for (std::size_t i = 0; i < a.inputs.size(); ++i) {
        expectSamePort(a.inputs[i], b.inputs[i]);
    }
    ASSERT_EQ(a.outputs.size(), b.outputs.size());
    for (std::size_t i = 0; i < a.outputs.size(); ++i) {
        expectSamePort(a.outputs[i], b.outputs[i]);
    }
    ASSERT_EQ(a.sections.size(), b.sections.size());
    for (std::size_t i = 0; i < a.sections.size(); ++i) {
        const Section& x = a.sections[i];
        const Section& y = b.sections[i];
        SCOPED_TRACE("section " + y.key);
        EXPECT_EQ(x.key, y.key);
        EXPECT_EQ(x.title, y.title);
        EXPECT_EQ(x.tooltip, y.tooltip);
        EXPECT_EQ(x.repeat, y.repeat);
        EXPECT_EQ(x.itemLabel, y.itemLabel);
        ASSERT_EQ(x.options.size(), y.options.size());
        for (std::size_t k = 0; k < x.options.size(); ++k) {
            expectSameOption(x.options[k], y.options[k]);
        }
        ASSERT_EQ(x.templates.size(), y.templates.size());
        for (std::size_t k = 0; k < x.templates.size(); ++k) {
            EXPECT_EQ(x.templates[k].id, y.templates[k].id);
            EXPECT_EQ(x.templates[k].label, y.templates[k].label);
            EXPECT_TRUE(same(Value(x.templates[k].values), Value(y.templates[k].values)))
                << y.templates[k].id;
        }
        EXPECT_EQ(x.defaultItems, y.defaultItems);
        EXPECT_EQ(x.requirements, y.requirements);
    }
}

std::string coerceError(const Option& option, const Value& input) {
    try {
        const Value out = option.coerce(input);
        return "no error, got " + pythonRepr(out);
    } catch (const std::invalid_argument& e) {
        return e.what();
    }
}

std::string message(const std::function<void()>& fn) {
    try {
        fn();
    } catch (const std::exception& e) {
        return e.what();
    }
    return "no error";
}

const std::size_t kBuiltinCount = 16;

}  // namespace

// ---- the registry (before any test registers further step types) -----------------------------

TEST(StepTypeRegistry, MatchesMdatkDump) {
    const Value dump = cases("step_types.json");
    const Array& tools = items(dump, "tools");
    ASSERT_EQ(tools.size(), kBuiltinCount);
    ASSERT_GE(stepTypes().size(), kBuiltinCount);
    for (std::size_t i = 0; i < tools.size(); ++i) {
        const StepType expected = stepTypeFrom(tools[i]);
        EXPECT_EQ(stepTypes()[i].name, expected.name);
        const StepType* found = findStepType(expected.name);
        ASSERT_NE(found, nullptr) << expected.name;
        EXPECT_EQ(found, &stepTypes()[i]);
        EXPECT_EQ(&stepType(expected.name), found);
        expectSameStepType(*found, expected);
    }
    EXPECT_EQ(stepCategories(), strings(*dump.find("categories")));
    const Array& kinds = items(dump, "file_kinds");
    ASSERT_EQ(fileKinds().size(), kinds.size());
    for (std::size_t i = 0; i < kinds.size(); ++i) {
        const Array& k = *kinds[i].get<Array>();
        EXPECT_EQ(fileKinds()[i].id, *k[0].get<std::string>());
        EXPECT_EQ(fileKinds()[i].label, *k[1].get<std::string>());
        EXPECT_EQ(fileKinds()[i].extensions, strings(k[2]));
    }
}

TEST(StepTypeRegistry, GetToolAndCategories) {
    ASSERT_EQ(stepTypes().size(), kBuiltinCount) << "run before the tests that register step types";
    const Value data = cases("registry_cases.json");
    for (const Value& c : items(data, "get_tool")) {
        const std::string name = str(c, "name");
        EXPECT_EQ(findStepType(name), nullptr);
        EXPECT_EQ(message([&] { stepType(name); }), str(c, "error"));
        try {
            stepType(name);
        } catch (const std::out_of_range&) {
            continue;
        }
        ADD_FAILURE() << "no std::out_of_range for " << name;
    }
    const auto groups = stepTypesByCategory();
    const Array& expected = items(data, "tools_by_category");
    ASSERT_EQ(groups.size(), expected.size());
    for (std::size_t i = 0; i < groups.size(); ++i) {
        const Array& e = *expected[i].get<Array>();
        EXPECT_EQ(groups[i].first, *e[0].get<std::string>());
        std::vector<std::string> names;
        for (const StepType* t : groups[i].second) {
            names.push_back(t->name);
        }
        EXPECT_EQ(names, strings(e[1]));
    }
}

TEST(StepTypeRegistry, ConcurrentFirstUse) {
    std::vector<std::thread> threads;
    std::vector<std::size_t> counts(8);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        threads.emplace_back(
            [&counts, i] { counts[i] = stepTypes().size() + (findStepType("vmd") ? 1 : 0); });
    }
    for (std::thread& t : threads) {
        t.join();
    }
    for (std::size_t n : counts) {
        EXPECT_EQ(n, stepTypes().size() + 1);
    }
}

// ---- schema ----------------------------------------------------------------------------------

TEST(StepTypeSchema, TemplatesAndDefaults) {
    const Value data = cases("registry_cases.json");
    for (const Value& c : items(data, "templates")) {
        const Section* section = stepType(str(c, "tool")).section(str(c, "section"));
        ASSERT_NE(section, nullptr);
        const std::string id = str(c, "id");
        if (const auto error = optionalString(*c.find("error"))) {
            EXPECT_EQ(message([&] { section->recordTemplate(id); }), *error);
            EXPECT_THROW(section->recordTemplate(id), std::out_of_range);
        } else {
            EXPECT_TRUE(same(Value(section->recordTemplate(id)), decode(*c.find("record")))) << id;
        }
    }
    const StepType& md = stepType("amber-md");
    EXPECT_EQ(md.section("nope"), nullptr);
    EXPECT_EQ(md.flatOption("name"), nullptr);  // an option of the repeated "stages" only
    ASSERT_NE(md.flatOption("engine"), nullptr);
    EXPECT_EQ(md.flatOption("engine")->label, "MD engine");
    EXPECT_EQ(md.section("stages")->option("name")->key, "name");
    EXPECT_EQ(md.section("stages")->option("engine"), nullptr);
    EXPECT_TRUE(same(*Value(md.defaults()).find("stages"),
                     Value(md.section("stages")->defaultRecords())));
}

TEST(StepTypeSchema, IssueTextAndKindNames) {
    EXPECT_EQ((Issue{"warning", "x: bad - using the default", {"x"}}).text(),
              "WARNING: x: bad - using the default");
    EXPECT_EQ((Issue{"error", "choose the file", {}}).text(), "ERROR: choose the file");
    EXPECT_EQ((Issue{"info", "", {}}).text(), "INFO: ");
    EXPECT_EQ(optionKindName(OptionKind::Multiline), "multiline");
    EXPECT_EQ(optionKindName(OptionKind::File), "file");
}

TEST(StepTypeSchema, KindsMatchMdatk) {
    const Value data = cases("registry_cases.json");
    for (const Value& c : items(data, "guess_kind")) {
        const Array& pair = *c.get<Array>();
        EXPECT_EQ(guessKind(*pair[0].get<std::string>()), *pair[1].get<std::string>())
            << pythonRepr(pair[0]);
    }
    for (const Value& c : items(data, "kind_label")) {
        const Array& pair = *c.get<Array>();
        EXPECT_EQ(kindLabel(*pair[0].get<std::string>()), *pair[1].get<std::string>());
    }
}

TEST(StepTypeSchema, RequirementsMatchMdatk) {
    const Value data = cases("requirements_cases.json");
    for (const Value& c : *data.get<Array>()) {
        const Value req = decode(*c.find("requires"));
        const Value values = decode(*c.find("values"));
        Requirements r(req.get<Table>()->begin(), req.get<Table>()->end());
        EXPECT_EQ(requirementsMet(r, *values.get<Table>()), *c.find("met")->get<bool>())
            << pythonRepr(req) << " / " << pythonRepr(values);
    }
}

// ---- Option::coerce --------------------------------------------------------------------------

namespace {

void expectCoerce(const Option& option, const Value& input, const Value& result) {
    SCOPED_TRACE(option.key + " <- " + pythonRepr(input));
    if (const Value* out = result.find("out")) {
        if (flag(result, "lossy")) {  // Python keeps a big integer: an error here
            EXPECT_EQ(coerceError(option, input),
                      option.key + ": expected an integer, got " + pythonRepr(input));
            return;
        }
        try {
            const Value got = option.coerce(input);
            EXPECT_TRUE(same(got, decode(*out)))
                << pythonRepr(got) << " != " << pythonRepr(decode(*out));
        } catch (const std::exception& e) {
            ADD_FAILURE() << "unexpected error " << e.what();
        }
        return;
    }
    std::string expected = str(result, "error");
    if (flag(result, "repr")) {
        expected += pythonRepr(input);
    }
    EXPECT_EQ(coerceError(option, input), expected);  // OverflowError (int(inf)) included
}

}  // namespace

TEST(StepTypeCoerce, MatchesMdatk) {
    const Value data = cases("coerce_cases.json");
    std::size_t count = 0;
    for (const Value& entry : items(data, "options")) {
        const Option option = optionFrom(*entry.find("option"));
        std::vector<Value> inputs;
        for (const Value& v : items(data, "inputs_" + str(entry, "inputs"))) {
            inputs.push_back(decode(v));
        }
        for (const Value& v : items(entry, "extra")) {
            inputs.push_back(decode(v));
        }
        const Array& results = items(entry, "results");
        ASSERT_EQ(inputs.size(), results.size());
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            expectCoerce(option, inputs[i], results[i]);
            ++count;
        }
    }
    EXPECT_GT(count, 10000u);
}

TEST(StepTypeCoerce, UnicodeDigitsAndWhitespace) {
    Option option;
    option.key = "n";
    option.kind = OptionKind::Int;
    const Value data = cases("digit_cases.json");
    for (const Value& c : *data.get<Array>()) {
        expectCoerce(option, decode(*c.find("in")), c);
    }
}

TEST(StepTypeCoerce, Messages) {
    const Option& nproc = *stepType("amber-md").flatOption("nproc");
    EXPECT_EQ(coerceError(nproc, Value("0")), "nproc: must be >= 1");
    EXPECT_EQ(coerceError(nproc, Value(1e9)), "nproc: must be <= 100000");
    EXPECT_EQ(coerceError(nproc, Value(std::nan(""))), "cannot convert float NaN to integer");
    EXPECT_EQ(coerceError(nproc, Value(-HUGE_VAL)), "cannot convert float infinity to integer");
    EXPECT_EQ(coerceError(nproc, Value(true)), "nproc: expected a number, got True");
    EXPECT_EQ(coerceError(nproc, Value()), "nproc: a value is required");
    const Option& engine = *stepType("amber-md").flatOption("engine");
    EXPECT_EQ(coerceError(engine, Value("x")),
              "engine: 'x' is not one of pmemd.cuda, pmemd, pmemd.MPI");
    EXPECT_EQ(engine.coerce(Value("pmemd")), Value("pmemd"));
}

// ---- normalise -------------------------------------------------------------------------------

TEST(StepTypeNormalise, MatchesMdatk) {
    const Value data = cases("normalise_cases.json");
    for (const Value& c : *data.get<Array>()) {
        const StepType& type = stepType(str(c, "tool"));
        SCOPED_TRACE(type.name + ": " + str(c, "label"));
        const Value settingsValue = decode(*c.find("settings"));
        const Table& settings = *settingsValue.get<Table>();
        const Value expected = decode(*c.find("result"));

        std::vector<Issue> issues;
        std::vector<std::string> lines;
        const auto warn = [&lines](const std::string& line) { lines.push_back(line); };
        const Table result = type.normalise(settings, &issues, warn);
        EXPECT_TRUE(same(Value(result), expected)) << pythonRepr(Value(result));
        const Array& expectedIssues = items(c, "issues");
        ASSERT_EQ(issues.size(), expectedIssues.size());
        for (std::size_t i = 0; i < issues.size(); ++i) {
            const Array& e = *expectedIssues[i].get<Array>();
            EXPECT_EQ(issues[i].level, *e[0].get<std::string>());
            EXPECT_EQ(issues[i].message, *e[1].get<std::string>());
            EXPECT_EQ(issues[i].keys, strings(e[2]));
            EXPECT_EQ(issues[i].text(), *items(c, "issue_texts")[i].get<std::string>());
        }
        EXPECT_EQ(lines, strings(*c.find("log_with_issues")));

        lines.clear();
        const Table without = type.normalise(settings, nullptr, warn);
        EXPECT_TRUE(same(Value(without), expected));
        EXPECT_EQ(lines, strings(*c.find("log_without_issues")));
    }
}

TEST(StepTypeNormalise, WarningsGoToTheLogWithoutCallback) {
    mdbridge::test::LogCapture log;
    const Table result = stepType("note").normalise(Table{{"text", 5}, {"b", 1}, {"a", 2}});
    EXPECT_EQ(result, (Table{{"text", "5"}}));
    EXPECT_EQ(log.lines(),
              std::vector<std::string>{"WARNING: note: ignoring unknown settings a, b"});
}

// ---- settings comments -----------------------------------------------------------------------

TEST(StepTypeComments, MatchesMdatkSettingsDoc) {
    const Value data = cases("setting_comments.json");
    for (const Value& t : items(data, "synthetic")) {
        if (!findStepType(str(t, "name"))) {
            registerStepType(stepTypeFrom(t));
        }
    }
    std::size_t count = 0;
    for (const Value& c : items(data, "cases")) {
        const std::string tool = str(c, "tool");
        const std::string key = str(c, "key");
        EXPECT_EQ(settingComment(tool, key), optionalString(*c.find("comment")))
            << tool << " / " << key;
        ++count;
    }
    EXPECT_GT(count, 400u);
}

TEST(StepTypeComments, OptionTextEdgeCases) {
    const Value data = cases("option_text_cases.json");
    StepType type;
    type.name = "parity-option-text";
    Section plain;
    plain.key = "plain";
    std::vector<std::pair<std::string, std::string>> expected;  // key, text
    for (const Value& c : *data.get<Array>()) {
        const std::string key = "k" + std::to_string(expected.size());
        if (flag(c, "section")) {
            Section s;
            s.key = key;
            s.repeat = true;
            s.tooltip = str(c, "tooltip");
            type.sections.push_back(std::move(s));
        } else {
            Option o;
            o.key = key;
            o.label = str(c, "label");
            o.tooltip = str(c, "tooltip");
            for (const Value& choice : items(c, "choices")) {
                const Array& pair = *choice.get<Array>();
                o.choices.push_back({*pair[0].get<std::string>(), *pair[1].get<std::string>()});
            }
            plain.options.push_back(std::move(o));
        }
        expected.emplace_back(key, str(c, "text"));
    }
    type.sections.insert(type.sections.begin(), std::move(plain));
    if (!findStepType(type.name)) {
        registerStepType(std::move(type));
    }
    for (const auto& [key, text] : expected) {
        EXPECT_EQ(settingComment("parity-option-text", key), text) << key;
    }
}

// ---- register (last: adds step types) --------------------------------------------------------

TEST(StepTypeRegister, ChecksLikeMdatk) {
    const Value data = cases("registry_cases.json");
    for (const Value& c : items(data, "register")) {
        StepType type = stepTypeFrom(*c.find("type"));
        const std::string name = type.name;
        if (const auto error = optionalString(*c.find("error"))) {
            EXPECT_EQ(message([&] { registerStepType(type); }), *error);
            EXPECT_THROW(registerStepType(type), std::invalid_argument);
        } else if (!findStepType(name)) {
            registerStepType(std::move(type));
            ASSERT_NE(findStepType(name), nullptr);
            EXPECT_EQ(stepTypes().back().name, name);
        }
    }
    const auto groups = stepTypesByCategory();
    EXPECT_EQ(groups.back().first, "Other");
    EXPECT_EQ(groups.back().second.front()->category, "Parity");
}
