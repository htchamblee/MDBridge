// The Lua host of the workflow generators (src/workflow/lua; docs/specs/m7-plan.md 4.4, 4.5): the
// sandbox, the memory and instruction limits and their issue texts, errors with file and line,
// core::Value <-> Lua round trips that keep the integer and float subtypes, and the same output
// bytes under two string-hash seeds.

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>
#include <lua.h>

#include "mdbridge/core/json.hpp"
#include "mdbridge/core/value.hpp"
#include "mdbridge/workflow/lua/convert.hpp"
#include "mdbridge/workflow/lua/host.hpp"

namespace {

using mdbridge::core::Array;
using mdbridge::core::DateTime;
using mdbridge::core::Table;
using mdbridge::core::Value;
using namespace mdbridge::workflow::lua;

std::vector<Value> eval(State& state, std::string_view source) {
    std::vector<Value> values;
    for (const Ref& ref : state.run(source)) {
        values.push_back(state.value(ref));
    }
    return values;
}

Value eval1(State& state, std::string_view source) {
    std::vector<Value> values = eval(state, source);
    return values.empty() ? Value() : values.front();
}

// Runs `source` and returns the Error it raises (fails the test if none).
Error errorOf(State& state, std::string_view source) {
    try {
        state.run(source);
    } catch (const Error& error) {
        return error;
    }
    ADD_FAILURE() << "no error from: " << source;
    return Error(Error::Kind::Runtime, "", 0, "");
}

Options withLimits(std::size_t memoryBytes, std::uint64_t instructions) {
    Options options;
    options.limits.memoryBytes = memoryBytes;
    options.limits.instructions = instructions;
    return options;
}

std::string str(const Value& value) {
    return value.stringOr("<not a string>");
}

// Every kind of value, nested, with keys out of sorted order.
Value everyKind() {
    Table inner;
    inner.emplace_back("zeta", Value(std::int64_t{1}));
    inner.emplace_back("alpha", Value(1.0));
    inner.emplace_back("mid", Value());
    DateTime when;
    when.kind = DateTime::Kind::OffsetDateTime;
    when.year = 2026;
    when.month = 11;
    when.day = 28;
    when.hour = 14;
    when.minute = 32;
    when.second = 5;
    when.microsecond = 250000;
    when.offsetMinutes = 120;
    Table top;
    top.emplace_back("ints", Value(Array{Value(0), Value(-7), Value(std::numeric_limits<std::int64_t>::max()),
                                         Value(std::numeric_limits<std::int64_t>::min())}));
    top.emplace_back("floats", Value(Array{Value(0.5), Value(1.0), Value(-2.0), Value(1e300), Value(5e-324),
                                           Value(std::numeric_limits<double>::infinity())}));
    top.emplace_back("strings", Value(Array{Value(""), Value(std::string("a\0b", 3)), Value("H\xc3\xa9lix \xe2\x80\x94 5'")}));
    top.emplace_back("flags", Value(Array{Value(true), Value(false)}));
    top.emplace_back("nulls", Value(Array{Value(), Value(1), Value()}));
    top.emplace_back("empty_array", Value(Array{}));
    top.emplace_back("empty_table", Value(Table{}));
    top.emplace_back("records", Value(Array{Value(inner), Value(Array{Value(Array{})})}));
    top.emplace_back("when", Value(when));
    top.emplace_back("b", Value("second key, sorted first"));
    return Value(top);
}

// --- sandbox -------------------------------------------------------------------------------------

TEST(LuaHostSandbox, OnlyTheSandboxLibrariesExist) {
    State state("steps/sandbox.lua");
    const std::vector<Value> types = eval(state, R"(
        return type(mdb), type(string), type(table), type(math), type(utf8),
               type(os), type(io), type(debug), type(package), type(coroutine), type(require),
               type(dofile), type(loadfile), type(print), type(warn), type(collectgarbage), type(rawset),
               type(string.dump), type(getmetatable("").__index.dump), type(loadstring), type(unpack)
    )");
    const std::vector<std::string> expected = {"table", "table", "table", "table", "table", "nil", "nil",
                                               "nil",   "nil",   "nil",   "nil",   "nil",   "nil", "nil",
                                               "nil",   "nil",   "nil",   "nil",   "nil",   "nil", "nil"};
    ASSERT_EQ(types.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(str(types[i]), expected[i]) << "result " << i + 1;
    }
}

TEST(LuaHostSandbox, BaseFunctionsAndLibrariesWork) {
    State state("steps/sandbox.lua");
    const std::vector<Value> values = eval(state, R"(
        local names = { "assert", "error", "ipairs", "pairs", "next", "pcall", "xpcall", "select",
                        "tonumber", "tostring", "type", "load", "getmetatable", "setmetatable",
                        "rawget", "rawequal", "rawlen" }
        local missing = {}
        for _, name in ipairs(names) do
            if type(_G[name]) ~= "function" then missing[#missing + 1] = name end
        end
        return table.concat(missing, ","), _G == _ENV, ("%d-%s"):format(7, "x"), utf8.len("h\u{e9}lix"),
               math.max(3, 9), table.unpack({ 4, 5 }), select("#", 1, 2, 3), _VERSION
    )");
    ASSERT_EQ(values.size(), 8u);  // table.unpack gives one value in the middle of the list
    EXPECT_EQ(str(values[0]), "");
    EXPECT_EQ(values[1], Value(true));
    EXPECT_EQ(str(values[2]), "7-x");
    EXPECT_EQ(values[3], Value(5));
    EXPECT_EQ(values[4], Value(9));
    EXPECT_EQ(values[5], Value(4));
    EXPECT_EQ(values[6], Value(3));
    EXPECT_EQ(str(values[7]), "Lua 5.4");
}

TEST(LuaHostSandbox, LoadReadsTextOnly) {
    State state("steps/load.lua");
    const std::vector<Value> values = eval(state, R"(
        local f = load("return 1 + 1")
        local g, message = load(string.char(27) .. "LuaT\0", "=binary")
        local h, message2 = load(string.char(27) .. "LuaT\0", "=binary", "b")
        local e = load("return x", "=env", "t", { x = 5 })
        local ok, bad = pcall(load, function() return nil end)
        return f(), g, message, h, message2, e(), ok, bad
    )");
    ASSERT_EQ(values.size(), 8u);
    EXPECT_EQ(values[0], Value(2));
    EXPECT_TRUE(values[1].isNull());
    EXPECT_NE(str(values[2]).find("attempt to load a binary chunk (mode is 't')"), std::string::npos) << str(values[2]);
    EXPECT_TRUE(values[3].isNull());  // asking for mode "b" does not help
    EXPECT_NE(str(values[4]).find("binary chunk"), std::string::npos);
    EXPECT_EQ(values[5], Value(5));
    EXPECT_EQ(values[6], Value(false));
    EXPECT_NE(str(values[7]).find("string expected (the sandbox loads Lua text only)"), std::string::npos)
        << str(values[7]);
}

TEST(LuaHostSandbox, MainChunkMustBeText) {
    State state("steps/binary.lua");
    const Error error = errorOf(state, std::string_view("\x1bLuaT\0\x19\x93\r\n", 10));
    EXPECT_EQ(error.kind(), Error::Kind::Syntax);
    EXPECT_EQ(error.file(), "steps/binary.lua");
    EXPECT_NE(error.message().find("attempt to load a binary chunk"), std::string::npos) << error.what();
}

TEST(LuaHostSandbox, RandomStartsFromAFixedSeed) {
    const char* source = R"(
        local t = {}
        for i = 1, 5 do t[i] = math.random(1, 1000000) end
        t[6] = math.random()
        return t
    )";
    State a("steps/random.lua");
    State b("steps/random.lua", Options{Limits{}, 12345u});  // the hash seed does not matter
    EXPECT_EQ(eval1(a, source), eval1(b, source));
    EXPECT_EQ(eval1(a, "math.randomseed(42) return math.random(1, 1000000)"),
              eval1(b, "math.randomseed(42) return math.random(1, 1000000)"));
    const Error error = errorOf(a, "math.randomseed()");
    EXPECT_EQ(error.kind(), Error::Kind::Runtime);
    EXPECT_EQ(error.line(), 1);
    EXPECT_EQ(error.message(), "math.randomseed needs a seed (the sandbox does not seed from the clock)");
}

// --- limits --------------------------------------------------------------------------------------

TEST(LuaHostLimits, AnEndlessLoopRunsOutOfTime) {
    State state("steps/loop.lua", withLimits(Limits{}.memoryBytes, 1'000'000));
    const Error error = errorOf(state, "local n = 0\nwhile true do\n  n = n + 1\nend\n");
    EXPECT_EQ(error.kind(), Error::Kind::OutOfTime);
    EXPECT_EQ(error.file(), "steps/loop.lua");
    EXPECT_TRUE(error.line() == 2 || error.line() == 3) << error.what();
    EXPECT_EQ(error.message(), "the generator ran out of time");
    EXPECT_EQ(error.issueText(), "internal error while checking this step: the generator ran out of time");
}

TEST(LuaHostLimits, TheDefaultBudgetEndsAnEndlessLoop) {
    State state("steps/loop.lua");
    EXPECT_EQ(errorOf(state, "while true do end").kind(), Error::Kind::OutOfTime);
}

TEST(LuaHostLimits, PcallAndXpcallCannotCatchRunningOutOfTime) {
    State state("steps/loop.lua", withLimits(Limits{}.memoryBytes, 1'000'000));
    EXPECT_EQ(errorOf(state, "local ok = pcall(function() while true do end end)\nreturn ok").kind(),
              Error::Kind::OutOfTime);
    EXPECT_EQ(errorOf(state, "for i = 1, 3 do pcall(function() while true do end end) end return 1").kind(),
              Error::Kind::OutOfTime);
    EXPECT_EQ(errorOf(state, "return xpcall(function() while true do end end, function(e) return e end)").kind(),
              Error::Kind::OutOfTime);
    EXPECT_EQ(errorOf(state, "return xpcall(error, function(e) while true do end end)").kind(), Error::Kind::OutOfTime);
    EXPECT_EQ(errorOf(state, "local t = setmetatable({}, { __index = function() while true do end end })\n"
                             "return pcall(function() return t.x end)")
                  .kind(),
              Error::Kind::OutOfTime);
}

TEST(LuaHostLimits, TheBudgetIsPerCall) {
    State state("steps/budget.lua", withLimits(Limits{}.memoryBytes, 1'000'000));
    const std::vector<Ref> results = state.run(R"(
        return function()
            local n = 0
            for i = 1, 100000 do n = n + i end  -- about 400 000 instructions
            return n
        end
    )");
    ASSERT_EQ(results.size(), 1u);
    for (int call = 0; call < 5; ++call) {  // 2 000 000 instructions in all
        const std::vector<Value> values = state.call(results[0], {});
        ASSERT_EQ(values.size(), 1u);
        EXPECT_EQ(values[0], Value(std::int64_t{5000050000}));
    }
}

TEST(LuaHostLimits, AGrowingTableRunsOutOfMemory) {
    State state("steps/memory.lua", withLimits(std::size_t{16} << 20, Limits{}.instructions));
    const Error error = errorOf(state, "local t = {}\nfor i = 1, 1e8 do t[i] = { i } end\n");
    EXPECT_EQ(error.kind(), Error::Kind::OutOfMemory);
    EXPECT_EQ(error.file(), "steps/memory.lua");
    EXPECT_EQ(error.message(), "the generator ran out of memory");
    EXPECT_EQ(error.issueText(), "internal error while checking this step: the generator ran out of memory");
    EXPECT_LE(state.memoryInUse(), std::size_t{16} << 20);
}

TEST(LuaHostLimits, OneLargeStringRunsOutOfMemory) {
    State state("steps/memory.lua");  // 256 MiB
    EXPECT_EQ(errorOf(state, "return string.rep('x', 1 << 30)").kind(), Error::Kind::OutOfMemory);
    EXPECT_EQ(errorOf(state, "return ('ab'):rep(1 << 28, ',')").kind(), Error::Kind::OutOfMemory);
}

TEST(LuaHostLimits, PcallAndLoadCannotCatchRunningOutOfMemory) {
    State state("steps/memory.lua", withLimits(std::size_t{8} << 20, Limits{}.instructions));
    EXPECT_EQ(errorOf(state, "local ok = pcall(string.rep, 'x', 1 << 30)\nreturn ok").kind(), Error::Kind::OutOfMemory);
    EXPECT_EQ(errorOf(state, "return xpcall(string.rep, function(e) return e end, 'x', 1 << 30)").kind(),
              Error::Kind::OutOfMemory);
    EXPECT_EQ(errorOf(state, "local source = 'return {' .. string.rep('1,', 1000000) .. '}'\n"
                             "return pcall(function() return load(source) end)")
                  .kind(),
              Error::Kind::OutOfMemory);
}

TEST(LuaHostLimits, TheStateCanBeUsedAfterALimit) {
    State state("steps/limits.lua", withLimits(std::size_t{16} << 20, 1'000'000));
    EXPECT_EQ(errorOf(state, "local t = {} for i = 1, 1e7 do t[i] = { i } end").kind(), Error::Kind::OutOfMemory);
    EXPECT_EQ(eval1(state, "return 1 + 1"), Value(2));
    EXPECT_EQ(errorOf(state, "while true do end").kind(), Error::Kind::OutOfTime);
    EXPECT_EQ(eval1(state, "local s = 0 for i = 1, 1000 do s = s + i end return s"), Value(500500));
}

TEST(LuaHostLimits, AStateNeedsMemoryForItsLibraries) {
    try {
        State state("steps/tiny.lua", withLimits(1024, Limits{}.instructions));
        FAIL() << "a state in 1 KiB";
    } catch (const Error& error) {
        EXPECT_EQ(error.kind(), Error::Kind::OutOfMemory);
    }
}

TEST(LuaHostLimits, FinalizersAreRefused) {
    // Lua runs __gc with hooks off, so the budget could not stop one that loops
    State state("steps/finalizer.lua", withLimits(Limits{}.memoryBytes, 1'000'000));
    const Error error = errorOf(state, "\nsetmetatable({}, { __gc = function() while true do end end })");
    EXPECT_EQ(error.kind(), Error::Kind::Runtime);
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "__gc metamethods are not available in the sandbox");
    // a __gc added later is never called (Lua marks objects for finalization in setmetatable)
    EXPECT_EQ(eval1(state, "local mt = {} local t = setmetatable({}, mt) mt.__gc = function() while true do end end\n"
                           "t = nil return getmetatable(setmetatable({}, { __index = {} })) ~= nil"),
              Value(true));
    // ~State closes the state: no finalizer runs
}

// --- errors --------------------------------------------------------------------------------------

TEST(LuaHostErrors, ARuntimeErrorCarriesFileAndLine) {
    State state("steps/errors.lua");
    const Error error = errorOf(state, "local x = 1\nlocal y = 2\nerror('boom')\n");
    EXPECT_EQ(error.kind(), Error::Kind::Runtime);
    EXPECT_EQ(error.file(), "steps/errors.lua");
    EXPECT_EQ(error.line(), 3);
    EXPECT_EQ(error.message(), "boom");
    EXPECT_STREQ(error.what(), "steps/errors.lua:3: boom");
    EXPECT_EQ(error.issueText(), "internal error while checking this step: steps/errors.lua:3: boom");
}

TEST(LuaHostErrors, LibraryAndVmErrorsCarryTheLine) {
    State state("steps/errors.lua");
    Error error = errorOf(state, "local s = 'a'\nreturn s:rep({})\n");
    EXPECT_EQ(error.line(), 2);
    EXPECT_NE(error.message().find("bad argument #1 to 'rep'"), std::string::npos) << error.what();
    error = errorOf(state, "\n\nreturn missing + 1\n");
    EXPECT_EQ(error.line(), 3);
    EXPECT_EQ(error.message(), "attempt to perform arithmetic on a nil value (global 'missing')");
    error = errorOf(state, "local t = nil\nreturn t.x\n");
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "attempt to index a nil value (local 't')");
}

TEST(LuaHostErrors, ErrorLevelTwoNamesTheCaller) {
    State state("steps/errors.lua");
    const Error error = errorOf(state, "local function need(x)\n  if not x then error('missing', 2) end\nend\n\nneed(nil)\n");
    EXPECT_EQ(error.line(), 5);
    EXPECT_EQ(error.message(), "missing");
}

TEST(LuaHostErrors, ARethrownErrorKeepsItsLine) {
    State state("steps/errors.lua");
    const Error error =
        errorOf(state, "local ok, e = pcall(function()\n  error('inner')\nend)\n\nerror(e, 0)\n");
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "inner");
}

TEST(LuaHostErrors, NonStringErrorObjects) {
    State state("steps/errors.lua");
    Error error = errorOf(state, "\nerror({ code = 1 })\n");
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "(error object is a table value)");
    error = errorOf(state, "error(setmetatable({}, { __tostring = function() return 'custom' end }))");
    EXPECT_EQ(error.line(), 1);
    EXPECT_EQ(error.message(), "custom");
}

TEST(LuaHostErrors, SyntaxErrors) {
    State state("steps/syntax.lua");
    const Error error = errorOf(state, "local a = 1\nlocal x = = 1\n");
    EXPECT_EQ(error.kind(), Error::Kind::Syntax);
    EXPECT_EQ(error.file(), "steps/syntax.lua");
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "unexpected symbol near '='");
    EXPECT_STREQ(error.what(), "steps/syntax.lua:2: unexpected symbol near '='");
}

TEST(LuaHostErrors, LongFileNamesAreKeptWhole) {
    const std::string file = "/home/someone/projects/a-very-long-folder-name/with/many/levels/steps/my_custom_step.lua";
    State state(file);
    Error error = errorOf(state, "\nerror('long')\n");
    EXPECT_EQ(error.file(), file);
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "long");
    error = errorOf(state, "\n\nlocal = 1\n");
    EXPECT_EQ(error.file(), file);
    EXPECT_EQ(error.line(), 3);
}

TEST(LuaHostErrors, AnErrorInAFunctionCalledFromCpp) {
    State state("steps/generator.lua");
    const std::vector<Ref> results = state.run(R"(local M = {}
function M.check(s)
  if s.count > 2 then
    error("too many: " .. s.count)
  end
  return "fine"
end
return M
)");
    ASSERT_EQ(results.size(), 1u);
    const Ref check = state.field(results[0], "check");
    ASSERT_TRUE(check);
    Table settings;
    settings.emplace_back("count", Value(1));
    EXPECT_EQ(state.call(check, std::vector<Value>{Value(settings)}), std::vector<Value>{Value("fine")});
    settings[0].second = Value(3);
    try {
        state.call(check, std::vector<Value>{Value(settings)});
        FAIL() << "no error";
    } catch (const Error& error) {
        EXPECT_EQ(error.line(), 4);
        EXPECT_EQ(error.message(), "too many: 3");
    }
}

TEST(LuaHostErrors, ErrorsOfLoadedChunksNameTheChunk) {
    State state("steps/errors.lua");
    const Error error = errorOf(state, "local f = load('\\n error(\"in chunk\")', '=helper')\nf()\n");
    EXPECT_EQ(error.file(), "helper");
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "in chunk");
}

int throwingFunction(lua_State* L) {
    return guarded(L, [&]() -> int { throw std::runtime_error("the C++ side failed"); });
}

int returningFunction(lua_State* L) {
    return guarded(L, [&] {
        lua_pushinteger(L, 41 + lua_tointeger(L, 1));
        return 1;
    });
}

TEST(LuaHostErrors, CppExceptionsInCFunctionsBecomeLuaErrors) {
    State state("steps/cfunctions.lua");
    state.protect([](lua_State* L) {
        lua_getglobal(L, "mdb");
        lua_pushcfunction(L, throwingFunction);
        lua_setfield(L, -2, "fail");
        lua_pushcfunction(L, returningFunction);
        lua_setfield(L, -2, "answer");
    });
    EXPECT_EQ(eval1(state, "return mdb.answer(1)"), Value(42));
    Error error = errorOf(state, "local x = 1\nmdb.fail()\n");
    EXPECT_EQ(error.kind(), Error::Kind::Runtime);
    EXPECT_EQ(error.line(), 2);
    EXPECT_EQ(error.message(), "the C++ side failed");
    const std::vector<Value> caught = eval(state, "return pcall(mdb.fail)");
    ASSERT_EQ(caught.size(), 2u);
    EXPECT_EQ(caught[0], Value(false));
    EXPECT_EQ(str(caught[1]), "the C++ side failed");  // called from pcall itself: no Lua line
}

TEST(LuaHostErrors, CppExceptionsInProtectPropagate) {
    State state("steps/protect.lua");
    EXPECT_THROW(state.protect([](lua_State*) { throw std::invalid_argument("bad"); }), std::invalid_argument);
    EXPECT_THROW(state.protect([&](lua_State*) { state.protect([](lua_State*) {}); }), std::logic_error);
    EXPECT_EQ(eval1(state, "return 'still fine'"), Value("still fine"));
}

// --- conversion ----------------------------------------------------------------------------------

TEST(LuaHostConversion, RoundTripsKeepEveryValue) {
    const Value value = everyKind();
    State state("steps/convert.lua");
    // mutable tables, straight through the C API
    Value back;
    state.protect([&](lua_State* L) {
        valueToLua(L, value, Access::Mutable);
        back = luaToValue(L, -1);
    });
    EXPECT_EQ(back, value);
    EXPECT_EQ(mdbridge::core::toJson(back), mdbridge::core::toJson(value));
    // read-only views, through a Lua function and back
    const std::vector<Ref> identity = state.run("return function(...) return ... end");
    const std::vector<Value> results = state.call(identity[0], std::vector<Value>{value, Value(Array{}), Value(Table{})});
    ASSERT_EQ(results.size(), 3u);
    EXPECT_EQ(results[0], value);
    EXPECT_EQ(results[1], Value(Array{}));
    EXPECT_EQ(results[2], Value(Table{}));
    // tables rebuilt in Lua from the view keep the order of pairs (the Value's)
    const std::vector<Ref> copy = state.run(R"(
        local function copy(x)
            if type(x) ~= "table" then return x end
            local t = {}
            for k, v in pairs(x) do t[k] = copy(v) end
            return t
        end
        return copy
    )");
    const Value copied = state.call(copy[0], std::vector<Value>{value})[0];
    EXPECT_EQ(*copied.find("ints"), *value.find("ints"));
    EXPECT_EQ(copied.find("records")->get<Array>()->at(0).get<Table>()->size(), 3u);
}

TEST(LuaHostConversion, SpecialFloats) {
    State state("steps/convert.lua");
    const std::vector<Ref> identity = state.run("return function(...) return ... end");
    const std::vector<Value> results = state.call(
        identity[0], std::vector<Value>{Value(-0.0), Value(std::numeric_limits<double>::quiet_NaN()),
                                        Value(-std::numeric_limits<double>::infinity())});
    ASSERT_EQ(results.size(), 3u);
    ASSERT_NE(results[0].get<double>(), nullptr);
    EXPECT_TRUE(std::signbit(*results[0].get<double>()));
    ASSERT_NE(results[1].get<double>(), nullptr);
    EXPECT_TRUE(std::isnan(*results[1].get<double>()));
    EXPECT_EQ(results[2], Value(-std::numeric_limits<double>::infinity()));
}

TEST(LuaHostConversion, IntegerAndFloatSubtypes) {
    State state("steps/convert.lua");
    const std::vector<Ref> types = state.run(R"(
        return function(t)
            local out = {}
            for i, v in ipairs(t) do out[i] = math.type(v) end
            return out
        end
    )");
    const Value seen = state.call(types[0], std::vector<Value>{Value(Array{Value(1), Value(1.0), Value(2.5),
                                                                          Value(std::numeric_limits<std::int64_t>::max()),
                                                                          Value(-0.0)})})[0];
    EXPECT_EQ(seen, Value(Array{Value("integer"), Value("float"), Value("float"), Value("integer"), Value("float")}));
    const std::vector<Value> values = eval(state, "return 3, 3.0, 7 // 2, 7 / 7, 2^53, math.maxinteger, -0.0, 1e100");
    const std::vector<Value> expected = {Value(3),   Value(3.0),
                                         Value(3),   Value(1.0),
                                         Value(9007199254740992.0), Value(std::numeric_limits<std::int64_t>::max()),
                                         Value(-0.0), Value(1e100)};
    EXPECT_EQ(values, expected);
    EXPECT_EQ(values[1].type(), Value::Type::Float);
    EXPECT_EQ(values[2].type(), Value::Type::Int);
}

TEST(LuaHostConversion, PlainTablesConvertWithoutHashOrder) {
    State state("steps/convert.lua");
    Table sorted;
    sorted.emplace_back("a", Value(2));
    sorted.emplace_back("b", Value(1));
    sorted.emplace_back("c", Value(Array{Value(3), Value(4)}));
    EXPECT_EQ(eval1(state, "return { b = 1, c = { 3, 4 }, a = 2 }"), Value(sorted));
    EXPECT_EQ(eval1(state, "return {}"), Value(Array{}));
    EXPECT_EQ(eval1(state, "return { [2] = 'b', [1] = 'a' }"), Value(Array{Value("a"), Value("b")}));
    EXPECT_EQ(eval1(state, "return { 1, mdb.null, 3 }"), Value(Array{Value(1), Value(), Value(3)}));
    EXPECT_TRUE(eval1(state, "return mdb.null").isNull());
    EXPECT_TRUE(eval1(state, "return nil").isNull());
}

TEST(LuaHostConversion, ValuesWithoutAValueAreRefused) {
    State state("steps/convert.lua");
    const auto refused = [&](std::string_view source) -> std::string {
        const std::vector<Ref> made = state.run(std::string("return function() ") + std::string(source) + " end");
        try {
            state.call(made[0], {});
        } catch (const Error& error) {
            EXPECT_EQ(error.kind(), Error::Kind::Conversion) << error.what();
            EXPECT_EQ(error.file(), "steps/convert.lua");
            EXPECT_EQ(error.line(), 0);
            return error.message();
        }
        ADD_FAILURE() << "converted: " << source;
        return {};
    };
    EXPECT_EQ(refused("return { files = { 'a', function() end } }"),
              "the result.files[2]: cannot convert a function to a value");
    EXPECT_EQ(refused("return 1, { x = { y = load } }"), "result 2.x.y: cannot convert a function to a value");
    EXPECT_EQ(refused("return { 1, 2, nil, 4 }"), "the result: a table whose integer keys are not 1..3 cannot be converted to a value");
    EXPECT_EQ(refused("return { 1, x = 2 }"), "the result: a table with both integer and string keys cannot be converted to a value");
    EXPECT_EQ(refused("return { [1.5] = 1 }"), "the result: a table with keys of type number cannot be converted to a value");
    EXPECT_EQ(refused("return { [true] = 1 }"), "the result: a table with keys of type boolean cannot be converted to a value");
    EXPECT_EQ(refused("local t = {} t.self = t return t"),
              "the result.self: a table that contains itself cannot be converted to a value");
    std::string deep = "the result";
    for (int level = 0; level < kMaxDepth; ++level) {
        deep += "[1]";
    }
    EXPECT_EQ(refused("local t = {} for i = 1, 150 do t = { t } end return t"),
              deep + ": the value is nested more than 100 levels deep");
    EXPECT_EQ(refused("return { ['a key'] = function() end }"),
              "the result[\"a key\"]: cannot convert a function to a value");
}

TEST(LuaHostConversion, TablesChangedInLuaKeepTheirKindAndOrder) {
    State state("steps/convert.lua");
    Table table;
    table.emplace_back("zeta", Value(1));
    table.emplace_back("alpha", Value(2));
    table.emplace_back("mid", Value(3));
    Value result;
    state.protect([&](lua_State* L) {
        valueToLua(L, Value(table), Access::Mutable);
        lua_setglobal(L, "t");
        valueToLua(L, Value(Array{}), Access::Mutable);
        lua_setglobal(L, "list");
        valueToLua(L, Value(Table{}), Access::Mutable);
        lua_setglobal(L, "dict");
    });
    result = eval1(state, R"(
        t.alpha = nil
        t.new2 = 5
        t.new1 = 4
        t.zeta = 10
        return { t = t, list = list, dict = dict }
    )");
    Table expected;
    expected.emplace_back("zeta", Value(10));
    expected.emplace_back("mid", Value(3));
    expected.emplace_back("new1", Value(4));
    expected.emplace_back("new2", Value(5));
    EXPECT_EQ(*result.find("t"), Value(expected));
    EXPECT_EQ(*result.find("list"), Value(Array{}));
    EXPECT_EQ(*result.find("dict"), Value(Table{}));
    // an array stays an array, so string keys added to it are refused
    try {
        eval(state, "list.x = 1 return list");
        FAIL() << "converted";
    } catch (const Error& error) {
        EXPECT_EQ(error.message(), "the value: an array with string keys cannot be converted to a value");
    }
    // markKind: an empty table made in Lua as a table
    Value marked;
    state.protect([&](lua_State* L) {
        lua_createtable(L, 0, 0);
        markKind(L, -1, Value::Type::Table);
        marked = luaToValue(L, -1);
    });
    EXPECT_EQ(marked, Value(Table{}));
}

TEST(LuaHostConversion, ReadOnlyViews) {
    Table record1;
    record1.emplace_back("b", Value(1));
    Table settings;
    settings.emplace_back("name", Value("x"));
    settings.emplace_back("list", Value(Array{Value(1), Value(2), Value(3)}));
    settings.emplace_back("records", Value(Array{Value(record1)}));
    settings.emplace_back("alpha", Value(false));
    const Value s(settings);
    State state("steps/view.lua");
    const std::vector<Ref> checks = state.run(R"(
        return function(s)
            local keys = {}
            for k in pairs(s) do keys[#keys + 1] = k end
            local sum = 0
            for _, v in ipairs(s.list) do sum = sum + v end
            local okSet, setError = pcall(function() s.name = "y" end)
            local okNested, nestedError = pcall(function() s.list[1] = 5 end)
            local okNew = pcall(function() s.records[1].c = 1 end)
            local okMeta = pcall(setmetatable, s, {})
            return { table.concat(keys, ","), #s.list, sum, table.concat(s.list, "-"), okSet, setError,
                     okNested, nestedError, okNew, getmetatable(s), okMeta, s.records[1].b, s.missing == nil }
        end
    )");
    const Value seen = state.call(checks[0], std::vector<Value>{s})[0];
    const Array& a = *seen.get<Array>();
    ASSERT_EQ(a.size(), 13u);
    EXPECT_EQ(str(a[0]), "name,list,records,alpha");  // the Value's order, not hash order
    EXPECT_EQ(a[1], Value(3));
    EXPECT_EQ(a[2], Value(6));
    EXPECT_EQ(str(a[3]), "1-2-3");
    EXPECT_EQ(a[4], Value(false));
    EXPECT_EQ(str(a[5]), "steps/view.lua:7: attempt to modify a read-only table");
    EXPECT_EQ(a[6], Value(false));
    EXPECT_EQ(str(a[7]), "steps/view.lua:8: attempt to modify a read-only table");
    EXPECT_EQ(a[8], Value(false));
    EXPECT_EQ(str(a[9]), "read-only table");
    EXPECT_EQ(a[10], Value(false));
    EXPECT_EQ(a[11], Value(1));
    EXPECT_EQ(a[12], Value(true));
    // and a view converts back to the value it shows
    const std::vector<Ref> identity = state.run("return function(x) return x end");
    EXPECT_EQ(state.call(identity[0], std::vector<Value>{s})[0], s);
}

TEST(LuaHostConversion, DateTimes) {
    State state("steps/datetime.lua");
    const Value value = everyKind();
    const Value& when = *value.find("when");
    const std::vector<Ref> f = state.run("return function(w) return tostring(w), w == w, type(w), w end");
    const std::vector<Value> results = state.call(f[0], std::vector<Value>{when});
    ASSERT_EQ(results.size(), 4u);
    EXPECT_EQ(str(results[0]), "2026-11-28T14:32:05.250000+02:00");
    EXPECT_EQ(results[1], Value(true));
    EXPECT_EQ(str(results[2]), "userdata");
    EXPECT_EQ(results[3], when);
}

TEST(LuaHostConversion, SortedPairs) {
    State state("steps/sorted.lua");
    const Value value = eval1(state, R"(
        local t = { b = 1, a = 2, [3] = "x", [1] = "y", [2.5] = "z", ["10"] = 0 }
        local keys = {}
        for k, v in mdb.sorted_pairs(t) do keys[#keys + 1] = math.type(k) or k end
        local ok, message = pcall(mdb.sorted_pairs, { [true] = 1 })
        return { table.concat(keys, ","), ok, message }
    )");
    const Array& a = *value.get<Array>();
    EXPECT_EQ(str(a[0]), "integer,float,integer,10,a,b");
    EXPECT_EQ(a[1], Value(false));
    EXPECT_EQ(str(a[2]), "sorted_pairs: keys of type boolean cannot be sorted");
}

// --- string-hash seeds ---------------------------------------------------------------------------

const char* const kHashOrderScript = R"(
    local t = {}
    for i = 1, 40 do t["key" .. i] = i end
    local order = {}
    for k in pairs(t) do order[#order + 1] = k end  -- hash order, on purpose
    return table.concat(order, ",")
)";

TEST(LuaHostHashSeed, TheSeedChangesTheOrderOfPairs) {
    // proves that the seed reaches Lua's string hashes (else the next test would show nothing)
    std::vector<std::string> orders;
    for (const unsigned int seed : {1u, 2u, 3u, 0x9e3779b9u, kDefaultHashSeed}) {
        State state("steps/seed.lua", Options{Limits{}, seed});
        orders.push_back(str(eval1(state, kHashOrderScript)));
    }
    bool differs = false;
    for (const std::string& order : orders) {
        differs = differs || order != orders.front();
    }
    EXPECT_TRUE(differs) << orders.front();
    // and the default seed is fixed: the same order in every state
    State a("steps/seed.lua");
    State b("steps/seed.lua");
    EXPECT_EQ(eval1(a, kHashOrderScript), eval1(b, kHashOrderScript));
}

TEST(LuaHostHashSeed, AGeneratorGivesTheSameBytesUnderTwoSeeds) {
    const char* generator = R"(
        local M = {}
        function M.generate(s, words)
            local lines = { "# " .. s.title }
            for key, value in pairs(s.options) do  -- a read-only view: the Value's order
                lines[#lines + 1] = key .. " = " .. value
            end
            local counts = {}
            for _, word in ipairs(words) do counts[word] = (counts[word] or 0) + 1 end
            for word, n in mdb.sorted_pairs(counts) do
                lines[#lines + 1] = word .. ": " .. n
            end
            return { text = table.concat(lines, "\n") .. "\n", counts = counts, files = { { name = "a.in" } } }
        end
        return M
    )";
    Table options;
    for (const char* key : {"nstlim", "dt", "cut", "ntb", "ntp", "igb", "temp0", "gamma_ln"}) {
        options.emplace_back(key, Value(std::string(key) + "-value"));
    }
    Table settings;
    settings.emplace_back("title", Value("production"));
    settings.emplace_back("options", Value(options));
    Array words;
    for (const char* word : {"helix", "loop", "stem", "helix", "bulge", "stem", "helix", "pseudoknot", "junction",
                             "loop", "hairpin", "kissing", "triplex", "stem"}) {
        words.emplace_back(word);
    }
    std::vector<std::string> outputs;
    for (const unsigned int seed : {1u, 0x9e3779b9u}) {
        State state("steps/amber_like.lua", Options{Limits{}, seed});
        const std::vector<Ref> module = state.run(generator);
        const Ref generate = state.field(module[0], "generate");
        const Value result = state.call(generate, std::vector<Value>{Value(settings), Value(words)})[0];
        outputs.push_back(mdbridge::core::toJson(result));
    }
    EXPECT_EQ(outputs[0], outputs[1]);
    EXPECT_NE(outputs[0].find("# production\\nnstlim = nstlim-value\\ndt = dt-value\\ncut"), std::string::npos)
        << outputs[0];
    EXPECT_NE(outputs[0].find("bulge: 1\\nhairpin: 1\\nhelix: 3\\njunction: 1\\nkissing: 1\\nloop: 2\\n"), std::string::npos)
        << outputs[0];
}

// --- references ----------------------------------------------------------------------------------

TEST(LuaHostRefs, GlobalsFieldsAndLifetime) {
    Ref kept;
    {
        State state("steps/refs.lua");
        std::vector<Ref> results = state.run("answer = 42\nreturn { name = 'x' }, nil, 7");
        ASSERT_EQ(results.size(), 3u);
        EXPECT_TRUE(results[0]);
        EXPECT_FALSE(results[1]);  // nil
        EXPECT_EQ(state.typeName(results[0]), "table");
        EXPECT_EQ(state.typeName(results[1]), "nil");
        EXPECT_EQ(state.value(state.global("answer")), Value(42));
        EXPECT_FALSE(state.global("nothing"));
        EXPECT_EQ(state.value(state.field(results[0], "name")), Value("x"));
        EXPECT_FALSE(state.field(results[0], "other"));
        EXPECT_THROW(state.field(results[2], "x"), Error);
        State other("steps/other.lua");
        EXPECT_THROW(other.value(results[0]), std::invalid_argument);
        kept = std::move(results[0]);
        EXPECT_TRUE(kept);
    }
    EXPECT_FALSE(kept);  // its state is gone
}

}  // namespace
