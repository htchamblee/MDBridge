// core::Value <-> Lua (convert.hpp).
//
// Kinds and key orders live in two weak-keyed registry tables (table -> 1 array / 2 table,
// table -> sequence of its keys), so the tables a generator sees are plain Lua tables. A read-only
// view is an empty proxy whose metatable reads from the data table; a third weak-keyed table maps
// each proxy to its data, which is what luaToValue converts.

#include "mdbridge/workflow/lua/convert.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <new>
#include <string>
#include <utility>
#include <vector>

#include <lauxlib.h>
#include <lua.h>

#include "host_internal.hpp"
#include "mdbridge/workflow/lua/host.hpp"

namespace mdbridge::workflow::lua {

namespace {

// Registry keys (their addresses) and the null sentinel; not const, so that no linker folds them
// into one address (MSVC's /OPT:ICF may merge identical read-only data)
char kKindsKey = 0;
char kOrdersKey = 0;
char kViewsKey = 0;
char kNullSentinel = 0;
constexpr const char* kDateTimeMeta = "mdbridge.datetime";
constexpr const char* kReadOnlyMessage = "attempt to modify a read-only table";

constexpr lua_Integer kArrayKind = 1;
constexpr lua_Integer kTableKind = 2;

void* sentinel() { return &kNullSentinel; }

Error conversionError(lua_State* L, std::string message) {
    return Error(Error::Kind::Conversion, detail::stateFile(L), 0, std::move(message));
}

void ensureStack(lua_State* L, int slots) {
    if (lua_checkstack(L, slots) == 0) {
        throw Error(Error::Kind::OutOfMemory, detail::stateFile(L), 0, std::string(kOutOfMemory));
    }
}

int toInt(std::size_t n) {
    return n > static_cast<std::size_t>(std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max()
                                                                       : static_cast<int>(n);
}

// Pushes registry[key] (a weak-keyed table).
void pushRegistryTable(lua_State* L, const char* key) {
    lua_rawgetp(L, LUA_REGISTRYINDEX, key);
}

// registry[key][table at `index`] = value at the top (popped)
void setMark(lua_State* L, const char* key, int index) {
    index = lua_absindex(L, index);
    pushRegistryTable(L, key);
    lua_pushvalue(L, index);
    lua_pushvalue(L, -3);
    lua_rawset(L, -3);
    lua_pop(L, 2);
}

// Pushes registry[key][table at `index`] (nil if none).
int getMark(lua_State* L, const char* key, int index) {
    index = lua_absindex(L, index);
    pushRegistryTable(L, key);
    lua_pushvalue(L, index);
    const int type = lua_rawget(L, -2);
    lua_remove(L, -2);
    return type;
}

lua_Integer kindOf(lua_State* L, int index) {
    getMark(L, &kKindsKey, index);
    const lua_Integer kind = lua_isinteger(L, -1) ? lua_tointeger(L, -1) : 0;
    lua_pop(L, 1);
    return kind;
}

// If the value at `index` is a read-only view, replaces it on the stack by its data table.
void resolveView(lua_State* L, int index) {
    index = lua_absindex(L, index);
    if (lua_type(L, index) != LUA_TTABLE) {
        return;
    }
    if (getMark(L, &kViewsKey, index) == LUA_TTABLE) {
        lua_replace(L, index);
    } else {
        lua_pop(L, 1);
    }
}

// --- read-only views -------------------------------------------------------------------------------

int readOnlyNewIndex(lua_State* L) {
    return luaL_error(L, "%s", kReadOnlyMessage);
}

int viewLength(lua_State* L) {
    lua_pushinteger(L, static_cast<lua_Integer>(lua_rawlen(L, lua_upvalueindex(1))));
    return 1;
}

// The iterator of a view's __pairs: upvalues data, keys (a sequence, or false for an array) and the
// position reached
int viewNext(lua_State* L) {
    const lua_Integer position = lua_tointeger(L, lua_upvalueindex(3)) + 1;
    lua_pushinteger(L, position);
    lua_replace(L, lua_upvalueindex(3));
    if (lua_istable(L, lua_upvalueindex(2))) {
        if (lua_rawgeti(L, lua_upvalueindex(2), position) == LUA_TNIL) {
            return 1;
        }
        lua_pushvalue(L, -1);
        lua_rawget(L, lua_upvalueindex(1));
        return 2;
    }
    if (lua_rawgeti(L, lua_upvalueindex(1), position) == LUA_TNIL) {
        return 1;
    }
    lua_pushinteger(L, position);
    lua_insert(L, -2);
    return 2;
}

// __pairs of a view (upvalue: the data table): the keys in the Value's order
int viewPairs(lua_State* L) {
    lua_pushvalue(L, lua_upvalueindex(1));
    if (getMark(L, &kOrdersKey, -1) != LUA_TTABLE) {
        lua_pop(L, 1);
        lua_pushboolean(L, 0);
    }
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, viewNext, 3);
    lua_pushnil(L);
    lua_pushnil(L);
    return 3;
}

// Replaces the data table at the top by a read-only view of it.
void wrapReadOnly(lua_State* L) {
    ensureStack(L, 6);
    const int data = lua_gettop(L);
    lua_createtable(L, 0, 0);  // the proxy
    lua_createtable(L, 0, 5);  // its metatable
    lua_pushvalue(L, data);
    lua_setfield(L, -2, "__index");
    lua_pushcfunction(L, readOnlyNewIndex);
    lua_setfield(L, -2, "__newindex");
    lua_pushvalue(L, data);
    lua_pushcclosure(L, viewLength, 1);
    lua_setfield(L, -2, "__len");
    lua_pushvalue(L, data);
    lua_pushcclosure(L, viewPairs, 1);
    lua_setfield(L, -2, "__pairs");
    lua_pushliteral(L, "read-only table");
    lua_setfield(L, -2, "__metatable");
    lua_setmetatable(L, -2);
    lua_pushvalue(L, data);
    setMark(L, &kViewsKey, -2);  // views[proxy] = data
    lua_replace(L, data);
}

// --- DateTime ------------------------------------------------------------------------------------

int dateTimeToString(lua_State* L) {
    const auto* dt = static_cast<const core::DateTime*>(luaL_checkudata(L, 1, kDateTimeMeta));
    return guarded(L, [&] {
        const std::string text = dt->isoformat();
        lua_pushlstring(L, text.data(), text.size());
        return 1;
    });
}

int dateTimeEquals(lua_State* L) {
    const auto* a = static_cast<const core::DateTime*>(luaL_testudata(L, 1, kDateTimeMeta));
    const auto* b = static_cast<const core::DateTime*>(luaL_testudata(L, 2, kDateTimeMeta));
    lua_pushboolean(L, a != nullptr && b != nullptr && *a == *b);
    return 1;
}

// --- Value to Lua --------------------------------------------------------------------------------

void push(lua_State* L, const core::Value& value, Access access, int depth) {
    ensureStack(L, 8);
    switch (value.type()) {
        case core::Value::Type::Null:
            pushNull(L);
            return;
        case core::Value::Type::Bool:
            lua_pushboolean(L, *value.get<bool>() ? 1 : 0);
            return;
        case core::Value::Type::Int:
            lua_pushinteger(L, static_cast<lua_Integer>(*value.get<std::int64_t>()));
            return;
        case core::Value::Type::Float:
            lua_pushnumber(L, static_cast<lua_Number>(*value.get<double>()));
            return;
        case core::Value::Type::String: {
            const std::string& text = *value.get<std::string>();
            lua_pushlstring(L, text.data(), text.size());
            return;
        }
        case core::Value::Type::DateTime: {
            void* block = lua_newuserdatauv(L, sizeof(core::DateTime), 0);
            new (block) core::DateTime(*value.get<core::DateTime>());  // trivially destructible
            luaL_setmetatable(L, kDateTimeMeta);
            return;
        }
        case core::Value::Type::Array:
        case core::Value::Type::Table:
            break;
    }
    if (depth >= kMaxDepth) {
        throw conversionError(L, "the value is nested more than " + std::to_string(kMaxDepth) + " levels deep");
    }
    if (const core::Array* array = value.get<core::Array>()) {
        lua_createtable(L, toInt(array->size()), 0);
        lua_Integer index = 0;
        for (const core::Value& element : *array) {
            push(L, element, access, depth + 1);
            lua_rawseti(L, -2, ++index);
        }
        lua_pushinteger(L, kArrayKind);
        setMark(L, &kKindsKey, -2);
    } else {
        const core::Table& table = *value.get<core::Table>();
        lua_createtable(L, 0, toInt(table.size()));
        lua_createtable(L, toInt(table.size()), 0);  // the keys in order
        lua_Integer index = 0;
        for (const auto& [key, element] : table) {
            lua_pushlstring(L, key.data(), key.size());
            if (lua_rawget(L, -3) == LUA_TNIL) {  // a repeated key keeps its first place
                lua_pushlstring(L, key.data(), key.size());
                lua_rawseti(L, -3, ++index);
            }
            lua_pop(L, 1);
            lua_pushlstring(L, key.data(), key.size());
            push(L, element, access, depth + 1);
            lua_rawset(L, -4);
        }
        setMark(L, &kOrdersKey, -2);  // pops the keys
        lua_pushinteger(L, kTableKind);
        setMark(L, &kKindsKey, -2);
    }
    if (access == Access::ReadOnly) {
        wrapReadOnly(L);
    }
}

// --- Lua to Value --------------------------------------------------------------------------------

bool isIdentifier(const std::string& key) {
    if (key.empty() || (key[0] >= '0' && key[0] <= '9')) {
        return false;
    }
    return std::all_of(key.begin(), key.end(), [](char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
    });
}

std::string keySegment(const std::string& key) {
    if (isIdentifier(key)) {
        return "." + key;
    }
    std::string segment = "[\"";
    for (const char c : key) {
        if (c == '"' || c == '\\') {
            segment += '\\';
        }
        segment += c;
    }
    return segment + "\"]";
}

class ToValue {
public:
    ToValue(lua_State* L, std::string_view what) : m_L(L), m_path(what) {}

    core::Value convert(int index) {
        index = lua_absindex(m_L, index);
        switch (lua_type(m_L, index)) {
            case LUA_TNONE:
            case LUA_TNIL:
                return {};
            case LUA_TBOOLEAN:
                return core::Value(lua_toboolean(m_L, index) != 0);
            case LUA_TNUMBER:
                if (lua_isinteger(m_L, index)) {
                    return core::Value(static_cast<std::int64_t>(lua_tointeger(m_L, index)));
                }
                return core::Value(static_cast<double>(lua_tonumber(m_L, index)));
            case LUA_TSTRING: {
                std::size_t size = 0;
                const char* text = lua_tolstring(m_L, index, &size);
                return core::Value(std::string(text, size));
            }
            case LUA_TLIGHTUSERDATA:
                if (lua_touserdata(m_L, index) == sentinel()) {
                    return {};
                }
                break;
            case LUA_TUSERDATA:
                if (const auto* dt = static_cast<const core::DateTime*>(luaL_testudata(m_L, index, kDateTimeMeta))) {
                    return core::Value(*dt);
                }
                break;
            case LUA_TTABLE:
                return convertTable(index);
            default:
                break;
        }
        throw fail(std::string("cannot convert a ") + luaL_typename(m_L, index) + " to a value");
    }

private:
    Error fail(const std::string& message) const {
        return conversionError(m_L, m_path + ": " + message);
    }

    struct Keys {
        std::vector<std::string> strings;
        lua_Integer count = 0;
        lua_Integer integers = 0;
        lua_Integer minimum = std::numeric_limits<lua_Integer>::max();
        lua_Integer maximum = std::numeric_limits<lua_Integer>::min();
        const char* other = nullptr;  // the type of a key that is neither
    };

    Keys keysOf(int data) {
        Keys keys;
        lua_pushnil(m_L);
        while (lua_next(m_L, data) != 0) {
            lua_pop(m_L, 1);
            ++keys.count;
            if (lua_type(m_L, -1) == LUA_TSTRING) {
                std::size_t size = 0;
                const char* text = lua_tolstring(m_L, -1, &size);  // a string already: next stays valid
                keys.strings.emplace_back(text, size);
            } else if (lua_isinteger(m_L, -1)) {
                const lua_Integer key = lua_tointeger(m_L, -1);
                ++keys.integers;
                keys.minimum = std::min(keys.minimum, key);
                keys.maximum = std::max(keys.maximum, key);
            } else if (keys.other == nullptr) {
                keys.other = luaL_typename(m_L, -1);
            }
        }
        return keys;
    }

    core::Value convertTable(int index) {
        ensureStack(m_L, 8);
        lua_pushvalue(m_L, index);
        resolveView(m_L, -1);
        const int data = lua_gettop(m_L);
        const void* identity = lua_topointer(m_L, data);
        if (std::find(m_active.begin(), m_active.end(), identity) != m_active.end()) {
            throw fail("a table that contains itself cannot be converted to a value");
        }
        if (static_cast<int>(m_active.size()) >= kMaxDepth) {
            throw fail("the value is nested more than " + std::to_string(kMaxDepth) + " levels deep");
        }
        m_active.push_back(identity);
        const lua_Integer kind = kindOf(m_L, data);
        Keys keys = keysOf(data);
        if (keys.other != nullptr) {
            throw fail(std::string("a table with keys of type ") + keys.other + " cannot be converted to a value");
        }
        const bool hasStrings = !keys.strings.empty();
        const bool hasIntegers = keys.integers > 0;
        if (hasStrings && hasIntegers) {
            throw fail("a table with both integer and string keys cannot be converted to a value");
        }
        core::Value result;
        if (kind == kArrayKind || (kind != kTableKind && !hasStrings)) {
            if (hasStrings) {
                throw fail("an array with string keys cannot be converted to a value");
            }
            if (hasIntegers && (keys.minimum != 1 || keys.maximum != keys.count)) {
                throw fail("a table whose integer keys are not 1.." + std::to_string(keys.count) +
                           " cannot be converted to a value");
            }
            core::Array array;
            array.reserve(static_cast<std::size_t>(keys.count));
            for (lua_Integer i = 1; i <= keys.count; ++i) {
                lua_rawgeti(m_L, data, i);
                m_path += "[" + std::to_string(i) + "]";
                array.push_back(convert(-1));
                m_path.resize(m_path.size() - std::to_string(i).size() - 2);
                lua_pop(m_L, 1);
            }
            result = core::Value(std::move(array));
        } else {
            if (hasIntegers) {
                throw fail("a table with integer keys cannot be converted to a value with named keys");
            }
            core::Table entries;
            entries.reserve(keys.strings.size());
            std::vector<std::string> remaining = std::move(keys.strings);
            std::sort(remaining.begin(), remaining.end());
            std::vector<bool> used(remaining.size(), false);
            const auto add = [&](const std::string& key) {
                lua_pushlstring(m_L, key.data(), key.size());
                lua_rawget(m_L, data);
                const std::size_t length = m_path.size();
                m_path += keySegment(key);
                entries.emplace_back(key, convert(-1));
                m_path.resize(length);
                lua_pop(m_L, 1);
            };
            if (getMark(m_L, &kOrdersKey, data) == LUA_TTABLE) {  // the recorded order first
                const int order = lua_gettop(m_L);
                for (lua_Integer i = 1; lua_rawgeti(m_L, order, i) == LUA_TSTRING; ++i) {
                    std::size_t size = 0;
                    const char* text = lua_tolstring(m_L, -1, &size);
                    std::string key(text, size);
                    lua_pop(m_L, 1);
                    const auto at = std::lower_bound(remaining.begin(), remaining.end(), key);
                    if (at != remaining.end() && *at == key && !used[static_cast<std::size_t>(at - remaining.begin())]) {
                        used[static_cast<std::size_t>(at - remaining.begin())] = true;
                        add(key);
                    }
                }
                lua_pop(m_L, 1);  // the nil that ended the order
            }
            lua_pop(m_L, 1);  // the order table or nil
            for (std::size_t i = 0; i < remaining.size(); ++i) {
                if (!used[i]) {
                    add(remaining[i]);
                }
            }
            result = core::Value(std::move(entries));
        }
        m_active.pop_back();
        lua_pop(m_L, 1);  // data
        return result;
    }

    lua_State* m_L;
    std::string m_path;
    std::vector<const void*> m_active;
};

// --- mdb.sorted_pairs ----------------------------------------------------------------------------

struct SortKey {
    int type = LUA_TNUMBER;  // LUA_TNUMBER before LUA_TSTRING
    bool isInteger = false;
    lua_Integer integer = 0;
    lua_Number number = 0;
    std::string text;
};

bool operator<(const SortKey& a, const SortKey& b) {
    if (a.type != b.type) {
        return a.type == LUA_TNUMBER;
    }
    if (a.type == LUA_TSTRING) {
        return a.text < b.text;
    }
    if (a.isInteger && b.isInteger) {
        return a.integer < b.integer;
    }
    const lua_Number x = a.isInteger ? static_cast<lua_Number>(a.integer) : a.number;
    const lua_Number y = b.isInteger ? static_cast<lua_Number>(b.integer) : b.number;
    if (x != y) {
        return x < y;
    }
    return a.isInteger && !b.isInteger;  // equal values: the integer first
}

// upvalues: data, sorted keys, position
int sortedNext(lua_State* L) {
    const lua_Integer position = lua_tointeger(L, lua_upvalueindex(3)) + 1;
    lua_pushinteger(L, position);
    lua_replace(L, lua_upvalueindex(3));
    if (lua_rawgeti(L, lua_upvalueindex(2), position) == LUA_TNIL) {
        return 1;
    }
    lua_pushvalue(L, -1);
    lua_rawget(L, lua_upvalueindex(1));
    return 2;
}

}  // namespace

void pushNull(lua_State* L) {
    lua_pushlightuserdata(L, sentinel());
}

bool isNull(lua_State* L, int index) {
    return lua_type(L, index) == LUA_TLIGHTUSERDATA && lua_touserdata(L, index) == sentinel();
}

void valueToLua(lua_State* L, const core::Value& value, Access access) {
    push(L, value, access, 0);
}

core::Value luaToValue(lua_State* L, int index, std::string_view what) {
    return ToValue(L, what).convert(index);
}

void markKind(lua_State* L, int index, core::Value::Type kind) {
    ensureStack(L, 5);
    index = lua_absindex(L, index);
    lua_pushvalue(L, index);
    resolveView(L, -1);
    lua_pushinteger(L, kind == core::Value::Type::Table ? kTableKind : kArrayKind);
    setMark(L, &kKindsKey, -2);
    lua_pop(L, 1);
}

int sortedPairs(lua_State* L) {
    luaL_checktype(L, 1, LUA_TTABLE);
    lua_settop(L, 1);
    resolveView(L, 1);
    std::vector<SortKey> keys;
    const char* bad = nullptr;
    guarded(L, [&] {  // std::bad_alloc becomes a Lua error
        lua_pushnil(L);
        while (lua_next(L, 1) != 0) {
            lua_pop(L, 1);
            SortKey key;
            key.type = lua_type(L, -1);
            if (key.type == LUA_TSTRING) {
                std::size_t size = 0;
                const char* text = lua_tolstring(L, -1, &size);
                key.text.assign(text, size);
            } else if (key.type == LUA_TNUMBER) {
                key.isInteger = lua_isinteger(L, -1) != 0;
                key.integer = lua_tointeger(L, -1);
                key.number = lua_tonumber(L, -1);
            } else {
                bad = luaL_typename(L, -1);
                lua_pop(L, 1);
                break;
            }
            keys.push_back(std::move(key));
        }
        if (bad == nullptr) {
            std::sort(keys.begin(), keys.end());
        }
        return 0;
    });
    if (bad != nullptr) {
        return luaL_error(L, "sorted_pairs: keys of type %s cannot be sorted", bad);
    }
    lua_createtable(L, toInt(keys.size()), 0);
    lua_Integer index = 0;
    for (const SortKey& key : keys) {
        if (key.type == LUA_TSTRING) {
            lua_pushlstring(L, key.text.data(), key.text.size());
        } else if (key.isInteger) {
            lua_pushinteger(L, key.integer);
        } else {
            lua_pushnumber(L, key.number);
        }
        lua_rawseti(L, -2, ++index);
    }
    lua_pushvalue(L, 1);  // data
    lua_insert(L, -2);
    lua_pushinteger(L, 0);
    lua_pushcclosure(L, sortedNext, 3);
    lua_pushnil(L);
    lua_pushnil(L);
    return 3;
}

namespace detail {

void installConversion(lua_State* L) {
    for (char* key : {&kKindsKey, &kOrdersKey, &kViewsKey}) {
        lua_createtable(L, 0, 0);
        lua_createtable(L, 0, 1);
        lua_pushliteral(L, "k");
        lua_setfield(L, -2, "__mode");
        lua_setmetatable(L, -2);
        lua_rawsetp(L, LUA_REGISTRYINDEX, key);
    }
    luaL_newmetatable(L, kDateTimeMeta);
    lua_pushcfunction(L, dateTimeToString);
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, dateTimeEquals);
    lua_setfield(L, -2, "__eq");
    lua_pushliteral(L, "datetime");
    lua_setfield(L, -2, "__metatable");
    lua_pop(L, 1);
}

}  // namespace detail

}  // namespace mdbridge::workflow::lua
