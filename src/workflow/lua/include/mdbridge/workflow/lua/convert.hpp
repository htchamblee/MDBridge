#pragma once
// core::Value <-> Lua (docs/specs/m7-plan.md 4.4). Call these inside State::protect (host.hpp).
//
// Value to Lua:
//   Null      -> mdb.null (a light userdata; nil cannot be stored in arrays and tables)
//   Bool      -> boolean
//   Int       -> integer, Float -> float (math.type tells them apart: 1 and 1.0 stay different)
//   String    -> string (bytes, embedded zeros included)
//   DateTime  -> a userdata whose tostring is DateTime::isoformat()
//   Array     -> a sequence 1..n
//   Table     -> a table with string keys
// Arrays and tables remember their kind and, for tables, the order of their keys (outside the
// table, in the registry), so that a value comes back unchanged and an empty array stays an array.
// With Access::ReadOnly every array and table is a read-only view: assigning to it raises "attempt
// to modify a read-only table", `#`, ipairs, table.concat and indexing work, and `pairs` visits the
// keys in the Value's order.
//
// Lua to Value:
//   nil, mdb.null -> Null; booleans; integers -> Int; floats -> Float; strings -> String;
//   the DateTime userdata -> DateTime; tables (and read-only views) -> Array or Table:
//   - a table made by valueToLua keeps its kind; keys it still has come in their recorded order,
//     keys added in Lua follow in sorted order (bytes);
//   - any other table with only the keys 1..n is an Array (the empty table too), one with only
//     string keys a Table with its keys sorted by bytes, so the result never depends on Lua's hash
//     order (or seed);
//   - anything else (functions, other userdata, tables with other keys, holes, cycles, more than
//     kMaxDepth levels) raises Error(Kind::Conversion) naming where in the value it is.

#include <string_view>

#include "mdbridge/core/value.hpp"

struct lua_State;

namespace mdbridge::workflow::lua {

enum class Access { Mutable, ReadOnly };

// Levels of nested arrays and tables either way.
inline constexpr int kMaxDepth = 100;

// Pushes `value` onto the stack.
void valueToLua(lua_State* L, const core::Value& value, Access access = Access::Mutable);

// The value at `index`. `what` names it in errors ("the result", "argument 2", ...): the message is
// e.g. "the result[2].files: cannot convert a function to a value". Throws Error(Kind::Conversion)
// with the file of the state and line 0.
core::Value luaToValue(lua_State* L, int index, std::string_view what = "the value");

// Pushes mdb.null.
void pushNull(lua_State* L);
bool isNull(lua_State* L, int index);

// Marks the table at `index` as an array or a table for luaToValue (an empty table is an array
// otherwise). For the mdb API: e.g. an empty dict that a generator returns.
void markKind(lua_State* L, int index, core::Value::Type kind);

// mdb.sorted_pairs(t), a lua_CFunction: the iterator triple of a generic `for` over the table (or
// read-only view) `t` in key order: numbers ascending, then strings by bytes; other key types raise
// an error.
int sortedPairs(lua_State* L);

namespace detail {
// Called once per state by State (host.cpp): the registry tables and metatables conversion uses.
void installConversion(lua_State* L);
}  // namespace detail

}  // namespace mdbridge::workflow::lua
