// The determinism lint of the generators' Lua files (src/workflow/lua/lint.cpp; m7-plan 4.4).

#include <string>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/workflow/lua/lint.hpp"

namespace {

using namespace mdbridge::workflow::lua;

// "line:column rule" of each finding
std::vector<std::string> findings(std::string_view source) {
    std::vector<std::string> out;
    for (const LintFinding& finding : lint(source)) {
        out.push_back(std::to_string(finding.line) + ":" + std::to_string(finding.column) + " " + finding.rule);
    }
    return out;
}

using Findings = std::vector<std::string>;

TEST(LuaLint, ACleanGeneratorPasses) {
    const char* source = R"(-- a generator in the style the rules ask for
local M = {}
function M.generate(s, ctx)
    local lines = {}
    for i, stage in ipairs(s.stages) do
        lines[#lines + 1] = "stage " .. i .. ": " .. mdb.text.format("g", stage.dt)
    end
    for name, value in mdb.sorted_pairs(s.extra) do
        lines[#lines + 1] = name .. " = " .. mdb.text.repr(value)
    end
    local t = { pairs = 1, next = 2, print = 3 }  -- keys named like the functions
    local u = t.pairs + t.next + ctx:print() + #string.format("%d items, %s", 3, "x")
    local text = [[
pairs(x) and tostring(y) inside a long string ]] .. ("%5d|%-3s|%%"):format(1, "a")
    return text .. u .. 1 .. 0x1F
end
return M
)";
    EXPECT_EQ(findings(source), Findings{});
}

TEST(LuaLint, LibrariesOutsideTheSandbox) {
    EXPECT_EQ(findings("local f = io.open('x')\nos.time()\nprint(debug.traceback())\n"),
              (Findings{"1:11 library", "2:1 library", "3:1 library", "3:7 library"}));
    EXPECT_EQ(findings("require('x') dofile('a') loadfile('b') collectgarbage() rawset(t, 1, 2) package.path = ''"),
              (Findings{"1:1 library", "1:14 library", "1:26 library", "1:40 library", "1:57 library", "1:73 library"}));
    EXPECT_EQ(findings("coroutine.wrap(f) unpack(t) setfenv(1, {}) warn('x')"),
              (Findings{"1:1 library", "1:19 library", "1:29 library", "1:44 library"}));
    EXPECT_EQ(findings("x.os = 1 y:print() local t = { io = 1, debug = 2 } goto print ::os::"), Findings{});
}

TEST(LuaLint, PairsAndNext) {
    EXPECT_EQ(findings("for k, v in pairs(t) do end\nlocal k = next(t)\n"), (Findings{"1:13 pairs", "2:11 pairs"}));
    EXPECT_EQ(findings("for i, v in ipairs(t) do end for k in mdb.sorted_pairs(t) do end"), Findings{});
}

TEST(LuaLint, TostringAndFormats) {
    EXPECT_EQ(findings("local a = tostring(x)"), Findings{"1:11 tostring"});
    EXPECT_EQ(findings("string.format('%.3f', x)"), Findings{"1:1 format"});
    EXPECT_EQ(findings("string.format('%5.2e|%G|%a', x)"), Findings{"1:1 format"});
    EXPECT_EQ(findings("string.format(\"%d %s %x %q %% %5d %-3s\", 1, 'a', 2, 'b', 3, 'c')"), Findings{});
    EXPECT_EQ(findings("string.format 'n = %g'"), Findings{"1:1 format"});
    EXPECT_EQ(findings("string.format(fmt, x)"), Findings{"1:1 format"});
    EXPECT_EQ(findings("local f = string.format"), Findings{"1:11 format"});
    EXPECT_EQ(findings("local s = ('%.1f'):format(x)"), Findings{"1:20 format"});
    EXPECT_EQ(findings("local s = ('%d'):format(x)"), Findings{});
    EXPECT_EQ(findings("local s = fmt:format(x)"), Findings{"1:15 format"});
    EXPECT_EQ(findings("local s = mdb.text.format('.3f', x)"), Findings{});
}

TEST(LuaLint, ConcatenatedFloats) {
    EXPECT_EQ(findings("local s = 'dt = ' .. 0.002"), Findings{"1:19 concat"});
    EXPECT_EQ(findings("local s = 1e-3 .. ''"), Findings{"1:16 concat"});
    EXPECT_EQ(findings("local s = 'x' .. -2.5 .. 0x1p4"), (Findings{"1:15 concat", "1:23 concat"}));
    EXPECT_EQ(findings("local s = 'pi ' .. math.pi .. math.sqrt(2)"), (Findings{"1:17 concat", "1:28 concat"}));
    EXPECT_EQ(findings("local s = 'n ' .. 3 .. math.floor(x) .. 0x10 .. #t"), Findings{});
}

TEST(LuaLint, TheCollector) {
    EXPECT_EQ(findings("setmetatable(t, { __gc = f })\nlocal w = setmetatable({}, { __mode = 'k' })\nmt['__gc'] = f"),
              (Findings{"1:19 collector", "2:30 collector", "3:4 collector"}));
}

TEST(LuaLint, StringsAndCommentsAreNotCode) {
    EXPECT_EQ(findings("-- pairs(t) tostring(x) os.time()\nlocal s = \"pairs(t) .. 1.5\" .. 'io.open'\n"
                       "--[==[ print(x)\n os.exit() ]==]\nlocal l = [=[ pairs ]=]"),
              Findings{});
    EXPECT_EQ(findings("local s = 'a\\'b' .. \"c\\\"d\\\\\" .. 'e\\z\n   f'\nprint(s)"), Findings{"3:1 library"});
}

TEST(LuaLint, Waivers) {
    EXPECT_EQ(findings("for k in pairs(t) do end -- lint-allow: pairs (only counts keys)"), Findings{});
    EXPECT_EQ(findings("local a = tostring(pairs) -- lint-allow: pairs, tostring"), Findings{});
    EXPECT_EQ(findings("local a = tostring(pairs) -- lint-allow: pairs"), Findings{"1:11 tostring"});
    EXPECT_EQ(findings("-- lint-allow: pairs\nfor k in pairs(t) do end"), Findings{"2:10 pairs"});  // its own line only
    EXPECT_EQ(findings("local x = 1 -- lint-allow: sorting"), Findings{"1:13 waiver"});
}

TEST(LuaLint, UnfinishedStringsAndComments) {
    EXPECT_EQ(findings("local s = 'open\nprint(1)"), Findings{"1:11 syntax"});
    EXPECT_EQ(findings("pairs(t)\nlocal s = [[ never closed"), (Findings{"1:1 pairs", "2:11 syntax"}));
    EXPECT_EQ(findings("--[[ never closed"), Findings{"1:1 syntax"});
}

TEST(LuaLint, LineEndingsAndFormatting) {
    EXPECT_EQ(findings("local a = 1\r\nlocal b = tostring(a)\r\n"), Findings{"2:11 tostring"});
    const std::vector<LintFinding> found = lint("\n  print(1)");
    ASSERT_EQ(found.size(), 1u);
    EXPECT_EQ(formatFinding("steps/x.lua", found[0]),
              "steps/x.lua:2:3: library: 'print' does not exist in the sandbox (no os, io, debug, package or "
              "coroutine library, no files, printing or modules)");
}

}  // namespace
