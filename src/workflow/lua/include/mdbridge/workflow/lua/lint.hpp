#pragma once
// The determinism lint of the generators' Lua files (docs/specs/m7-plan.md 4.4, risk 1): byte-exact
// output must not depend on Lua's hash order, Lua's number formatting or anything outside the
// sandbox. It reads the tokens of a file (no Lua state) and reports:
//
//   library    os, io, debug, package, coroutine, require, dofile, loadfile, print, warn,
//              collectgarbage, rawset (and Lua 5.1's loadstring, setfenv, getfenv, module, unpack):
//              not in the sandbox at all
//   pairs      pairs and next: the order of a map's keys follows its hashes; use ipairs over arrays
//              or mdb.sorted_pairs
//   tostring   tostring: Lua writes floats unlike Python (%.14g); numbers go through mdb.text
//   format     string.format / :format with a float conversion (%a %e %f %g, either case) in its
//              literal format, or with a format that is not a literal string (it cannot be checked)
//   concat     `..` beside a float literal or a math function that gives floats (math.pi, sqrt ...)
//   collector  __gc and __mode: finalizers and weak tables follow the garbage collector's timing
//              (setmetatable refuses __gc at run time as well)
//   waiver     a waiver names an unknown rule
//   syntax     an unfinished string or comment (the rest of the file is not checked)
//
// Names are checked where they stand free (not after `.` or `:` and not as a key of a table
// constructor), so `x.print` and `{ next = 1 }` pass. A finding is waived by a comment on its line:
// `-- lint-allow: pairs` (several rules separated by spaces or commas), for a use whose order or
// formatting cannot reach the output. The lint is a guard, not a proof: a float that reaches `..` or
// string.format through a variable is not seen; every generator case also runs under two hash
// seeds (host.hpp).

#include <string>
#include <string_view>
#include <vector>

namespace mdbridge::workflow::lua {

struct LintFinding {
    int line = 0;    // 1-based
    int column = 0;  // 1-based, in bytes
    std::string rule;
    std::string message;
};

// The findings of `source` in the order of the file.
std::vector<LintFinding> lint(std::string_view source);

// "file:line:column: rule: message"
std::string formatFinding(std::string_view file, const LintFinding& finding);

}  // namespace mdbridge::workflow::lua
