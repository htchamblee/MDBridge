// The sandboxed Lua state of a generator file (host.hpp).

#include "mdbridge/workflow/lua/host.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>

#include "host_internal.hpp"
#include "lua_seed.hpp"
#include "mdbridge/workflow/lua/convert.hpp"

namespace mdbridge::workflow::lua {

namespace {

// What the allocator, the count hook and the sandbox functions share; the allocator's userdata (its
// SeedSource base gives lstate.c the hash seed, lua_seed.hpp).
struct Control : detail::SeedSource {
    std::string file;           // State::file()
    std::string chunkShortSrc;  // how Lua names the file in messages (lua_Debug::short_src)
    std::size_t used = 0;       // bytes allocated
    std::size_t limit = 0;      // Limits::memoryBytes
    std::uint64_t executed = 0;  // instructions of the current call
    std::uint64_t budget = 0;    // Limits::instructions
    bool active = false;         // inside State::protect
    bool closing = false;        // lua_close: finalizers get no budget
    bool timeExceeded = false;
    bool outOfMemory = false;    // a memory error seen (and passed on) by pcall / xpcall / load
    // where the budget ran out, and the error the message handler saw
    std::string timeFile;
    int timeLine = 0;
    bool handled = false;
    std::string errorFile;
    int errorLine = 0;
    std::string errorMessage;
};

Control& control(lua_State* L) {
    void* data = nullptr;
    lua_getallocf(L, &data);
    return *static_cast<Control*>(static_cast<detail::SeedSource*>(data));
}

void* allocate(void* data, void* block, std::size_t oldSize, std::size_t newSize) noexcept {
    Control& c = *static_cast<Control*>(static_cast<detail::SeedSource*>(data));
    const std::size_t old = block != nullptr ? oldSize : 0;  // without a block, oldSize is a type code
    if (newSize == 0) {
        std::free(block);
        c.used -= old;
        return nullptr;
    }
    if (newSize > old && newSize - old > c.limit - std::min(c.used, c.limit)) {
        return nullptr;  // over the cap: Lua collects garbage and tries again, then raises LUA_ERRMEM
    }
    void* grown = std::realloc(block, newSize);
    if (grown == nullptr) {
        return newSize <= old ? block : nullptr;  // Lua assumes that shrinking never fails
    }
    c.used = c.used - old + newSize;
    return grown;
}

// Where Lua names a function's chunk: the full name for files ("@...") and literal names ("=..."),
// Lua's short form for chunks loaded from strings.
std::string sourceName(const lua_Debug& ar) {
    if (ar.source != nullptr && (ar.source[0] == '@' || ar.source[0] == '=')) {
        return std::string(ar.source + 1);
    }
    return std::string(ar.short_src);
}

// If `message` starts with "<shortSrc>:<line>: ", that line and the rest of the message.
bool splitLocation(std::string_view message, std::string_view shortSrc, int& line, std::string& rest) {
    if (shortSrc.empty() || message.size() <= shortSrc.size() + 1 || message.substr(0, shortSrc.size()) != shortSrc ||
        message[shortSrc.size()] != ':') {
        return false;
    }
    std::size_t at = shortSrc.size() + 1;
    long long number = 0;
    const std::size_t digits = at;
    while (at < message.size() && message[at] >= '0' && message[at] <= '9' && number < 100'000'000) {
        number = number * 10 + (message[at] - '0');
        ++at;
    }
    if (at == digits || message.substr(at, 2) != ": ") {
        return false;
    }
    line = static_cast<int>(number);
    rest = std::string(message.substr(at + 2));
    return true;
}

void countHook(lua_State* L, lua_Debug* ar) {
    Control& c = control(L);
    c.executed += kInstructionStep;
    if (!c.closing && c.executed <= c.budget) {
        return;
    }
    if (!c.timeExceeded) {
        c.timeExceeded = true;
        try {
            if (lua_getinfo(L, "Sl", ar) != 0 && ar->currentline > 0) {
                c.timeFile = sourceName(*ar);
                c.timeLine = ar->currentline;
            }
        } catch (const std::exception&) {  // std::bad_alloc: no location
        }
    }
    lua_pushlstring(L, kOutOfTime.data(), kOutOfTime.size());
    lua_error(L);
}

// The message handler of State::protect: records the file, line and text of the error.
int messageHandler(lua_State* L) {
    Control& c = control(L);
    try {
        std::string message;
        if (lua_type(L, 1) == LUA_TSTRING || lua_type(L, 1) == LUA_TNUMBER) {
            std::size_t size = 0;
            lua_pushvalue(L, 1);
            const char* text = lua_tolstring(L, -1, &size);
            message.assign(text, size);
            lua_pop(L, 1);
        } else if (luaL_callmeta(L, 1, "__tostring") != 0 && lua_type(L, -1) == LUA_TSTRING) {
            message = lua_tostring(L, -1);
            lua_pop(L, 1);
        } else {
            message = std::string("(error object is a ") + luaL_typename(L, 1) + " value)";
        }
        std::string file = c.file;
        int line = 0;
        bool innermost = true;
        bool located = false;
        lua_Debug ar;
        for (int level = 1; !located && lua_getstack(L, level, &ar) != 0; ++level) {
            if (lua_getinfo(L, "Sl", &ar) == 0 || ar.currentline <= 0) {
                continue;
            }
            if (innermost) {  // by default the innermost Lua function
                file = sourceName(ar);
                line = ar.currentline;
                innermost = false;
            }
            std::string rest;
            int at = 0;
            if (splitLocation(message, ar.short_src, at, rest)) {  // "where:line: " of error(msg, level)
                file = sourceName(ar);
                line = at;
                message = std::move(rest);
                located = true;
            }
        }
        std::string rest;
        int at = 0;
        if (!located && splitLocation(message, c.chunkShortSrc, at, rest)) {  // a function that has returned
            file = c.file;
            line = at;
            message = std::move(rest);
        }
        c.errorFile = std::move(file);
        c.errorLine = line;
        c.errorMessage = std::move(message);
        c.handled = true;
    } catch (const std::exception&) {  // std::bad_alloc: the error is reported without its location
    }
    lua_settop(L, 1);
    return 1;
}

int panic(lua_State* L) {
    const char* message = lua_tostring(L, -1);
    std::fprintf(stderr, "MDBridge: Lua error outside a protected call: %s\n", message != nullptr ? message : "?");
    std::fflush(stderr);
    return 0;  // Lua aborts
}

// --- the sandbox's replacements -------------------------------------------------------------------

// After a protected call inside the sandbox: running out of time or memory is passed on.
void passLimits(lua_State* L, int status) {
    Control& c = control(L);
    if (status == LUA_ERRMEM) {
        c.outOfMemory = true;
    }
    if (status != LUA_OK && (c.outOfMemory || c.timeExceeded)) {
        lua_error(L);
    }
}

int sandboxPcall(lua_State* L) {
    luaL_checkany(L, 1);
    const int status = lua_pcall(L, lua_gettop(L) - 1, LUA_MULTRET, 0);
    passLimits(L, status);
    luaL_checkstack(L, 1, nullptr);  // the results may fill the stack
    lua_pushboolean(L, status == LUA_OK ? 1 : 0);
    lua_insert(L, 1);
    return lua_gettop(L);
}

// The message handler xpcall installs (upvalue: the caller's handler). Once the time or memory has
// run out, the caller's handler is not called: Lua runs a message handler before it unwinds, and an
// error raised by the count hook leaves hooks off until then, so a handler that loops could not be
// stopped.
int xpcallHandler(lua_State* L) {
    const Control& c = control(L);
    if (c.timeExceeded || c.outOfMemory) {
        lua_settop(L, 1);
        return 1;
    }
    lua_pushvalue(L, lua_upvalueindex(1));
    lua_insert(L, 1);
    lua_call(L, lua_gettop(L) - 1, 1);
    return 1;
}

int sandboxXpcall(lua_State* L) {
    const int n = lua_gettop(L);
    luaL_checktype(L, 2, LUA_TFUNCTION);
    lua_pushvalue(L, 2);
    lua_pushcclosure(L, xpcallHandler, 1);
    lua_replace(L, 2);
    lua_pushvalue(L, 1);  // f, h, args... -> h, f, args...
    lua_pushvalue(L, 2);
    lua_replace(L, 1);
    lua_replace(L, 2);
    const int status = lua_pcall(L, n - 2, LUA_MULTRET, 1);
    passLimits(L, status);
    luaL_checkstack(L, 1, nullptr);
    lua_pushboolean(L, status == LUA_OK ? 1 : 0);
    lua_replace(L, 1);
    return lua_gettop(L);
}

// load(chunk [, chunkname [, mode [, env]]]): strings of Lua text only (mode is always "t")
int sandboxLoad(lua_State* L) {
    if (lua_type(L, 1) != LUA_TSTRING) {
        return luaL_argerror(L, 1, "string expected (the sandbox loads Lua text only)");
    }
    std::size_t size = 0;
    const char* text = lua_tolstring(L, 1, &size);
    const char* chunkName = luaL_optstring(L, 2, text);
    (void)luaL_optstring(L, 3, "t");  // accepted for compatibility; the mode is always "t"
    const bool hasEnv = !lua_isnone(L, 4);
    const int status = luaL_loadbufferx(L, text, size, chunkName, "t");
    if (status == LUA_ERRMEM) {
        control(L).outOfMemory = true;
        return lua_error(L);
    }
    if (status != LUA_OK) {
        luaL_pushfail(L);
        lua_insert(L, -2);
        return 2;
    }
    if (hasEnv) {
        lua_pushvalue(L, 4);
        if (lua_setupvalue(L, -2, 1) == nullptr) {
            lua_pop(L, 1);
        }
    }
    return 1;
}

// setmetatable(t, mt) without finalizers (upvalue: Lua's setmetatable): Lua runs __gc metamethods
// with hooks off, outside the instruction budget, and at the collector's moments
int sandboxSetmetatable(lua_State* L) {
    if (lua_type(L, 2) == LUA_TTABLE) {
        lua_pushliteral(L, "__gc");
        if (lua_rawget(L, 2) != LUA_TNIL) {
            return luaL_error(L, "__gc metamethods are not available in the sandbox");
        }
        lua_pop(L, 1);
    }
    lua_pushvalue(L, lua_upvalueindex(1));
    lua_insert(L, 1);
    lua_call(L, lua_gettop(L) - 1, 1);
    return 1;
}

// math.randomseed(x [, y]) with an explicit seed only (upvalue: Lua's randomseed)
int sandboxRandomseed(lua_State* L) {
    if (lua_isnone(L, 1)) {
        return luaL_error(L, "math.randomseed needs a seed (the sandbox does not seed from the clock)");
    }
    lua_pushvalue(L, lua_upvalueindex(1));
    lua_insert(L, 1);
    lua_call(L, lua_gettop(L) - 1, LUA_MULTRET);
    return lua_gettop(L);
}

void openSandbox(lua_State* L, Control& c, const std::string& chunkName) {
    luaL_requiref(L, LUA_GNAME, luaopen_base, 1);
    luaL_requiref(L, LUA_STRLIBNAME, luaopen_string, 1);
    luaL_requiref(L, LUA_TABLIBNAME, luaopen_table, 1);
    luaL_requiref(L, LUA_MATHLIBNAME, luaopen_math, 1);
    luaL_requiref(L, LUA_UTF8LIBNAME, luaopen_utf8, 1);
    lua_settop(L, 0);
    for (const char* name : {"dofile", "loadfile", "print", "warn", "collectgarbage", "rawset"}) {
        lua_pushnil(L);
        lua_setglobal(L, name);
    }
    lua_pushcfunction(L, sandboxLoad);
    lua_setglobal(L, "load");
    lua_pushcfunction(L, sandboxPcall);
    lua_setglobal(L, "pcall");
    lua_pushcfunction(L, sandboxXpcall);
    lua_setglobal(L, "xpcall");
    lua_getglobal(L, "setmetatable");
    lua_pushcclosure(L, sandboxSetmetatable, 1);
    lua_setglobal(L, "setmetatable");

    lua_getglobal(L, LUA_STRLIBNAME);
    lua_pushnil(L);
    lua_setfield(L, -2, "dump");  // binary chunks cannot be loaded anyway
    lua_pop(L, 1);

    lua_getglobal(L, LUA_MATHLIBNAME);  // a fixed seed instead of the clock's and the state's address
    lua_getfield(L, -1, "randomseed");
    lua_pushvalue(L, -1);
    lua_pushinteger(L, 0);
    lua_call(L, 1, 0);
    lua_pushcclosure(L, sandboxRandomseed, 1);
    lua_setfield(L, -2, "randomseed");
    lua_pop(L, 1);

    detail::installConversion(L);
    lua_createtable(L, 0, 2);
    pushNull(L);
    lua_setfield(L, -2, "null");
    lua_pushcfunction(L, sortedPairs);
    lua_setfield(L, -2, "sorted_pairs");
    lua_setglobal(L, "mdb");

    // Lua's name for the file in messages (it shortens long names), for State::run's syntax errors
    if (luaL_loadbufferx(L, "return", 6, chunkName.c_str(), "t") == LUA_OK) {
        lua_Debug ar;
        lua_getinfo(L, ">S", &ar);
        c.chunkShortSrc = ar.short_src;
    }
    lua_settop(L, 0);
}

// The trampoline of State::protect: runs the body; a C++ exception is kept for the host (it must
// not unwind through Lua's own frames, which would turn it into an unknown error status)
struct ProtectedCall {
    const std::function<void(lua_State*)>* body = nullptr;
    std::exception_ptr exception;
};

int trampoline(lua_State* L) {
    auto* call = static_cast<ProtectedCall*>(lua_touserdata(L, 1));
    lua_settop(L, 0);
    try {
        (*call->body)(L);
    } catch (const std::exception&) {  // Lua's own errors are not std::exception and pass through
        call->exception = std::current_exception();
    }
    return 0;
}

int collectGarbage(lua_State* L) {
    lua_gc(L, LUA_GCCOLLECT);
    return 0;
}

std::string buildWhat(const std::string& file, int line, const std::string& message) {
    std::string text = file;
    if (line > 0) {
        text += ":" + std::to_string(line);
    }
    return text + ": " + message;
}

}  // namespace

// --- Error ---------------------------------------------------------------------------------------

Error::Error(Kind kind, std::string file, int line, std::string message)
    : std::runtime_error(buildWhat(file, line, message)),
      m_kind(kind),
      m_file(std::move(file)),
      m_line(line),
      m_message(std::move(message)) {}

std::string Error::issueText() const {
    switch (m_kind) {
        case Kind::OutOfTime:
            return std::string(kIssuePrefix) + std::string(kOutOfTime);
        case Kind::OutOfMemory:
            return std::string(kIssuePrefix) + std::string(kOutOfMemory);
        default:
            return std::string(kIssuePrefix) + what();
    }
}

// --- State ---------------------------------------------------------------------------------------

struct State::Impl {
    lua_State* L = nullptr;
    Control control;
    Options options;
    std::string chunkName;

    Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
    ~Impl() {
        if (L != nullptr) {
            control.closing = true;  // finalizers that loop are stopped by the count hook
            lua_close(L);
        }
    }
};

namespace detail {

const std::string& stateFile(lua_State* L) {
    return control(L).file;
}

}  // namespace detail

State::State(std::string fileName, Options options) : m_impl(std::make_shared<Impl>()) {
    Impl& impl = *m_impl;
    Control& c = impl.control;
    impl.options = options;
    c.file = std::move(fileName);
    c.hashSeed = options.hashSeed;
    c.limit = options.limits.memoryBytes;
    impl.chunkName = "@" + c.file;
    impl.L = lua_newstate(allocate, static_cast<detail::SeedSource*>(&c));
    if (impl.L == nullptr) {
        throw Error(Error::Kind::OutOfMemory, c.file, 0, std::string(kOutOfMemory));
    }
    lua_atpanic(impl.L, panic);
    lua_sethook(impl.L, countHook, LUA_MASKCOUNT, kInstructionStep);
    protect([&](lua_State* L) { openSandbox(L, c, impl.chunkName); });
}

State::~State() = default;

const std::string& State::file() const noexcept {
    return m_impl->control.file;
}

const Options& State::options() const noexcept {
    return m_impl->options;
}

std::size_t State::memoryInUse() const noexcept {
    return m_impl->control.used;
}

void State::protect(const std::function<void(lua_State*)>& body) {
    Impl& impl = *m_impl;
    Control& c = impl.control;
    lua_State* L = impl.L;
    if (c.active) {
        throw std::logic_error("lua::State::protect: already inside a protected call of this state");
    }
    struct Active {
        Control& c;
        ~Active() { c.active = false; }
    } active{c};
    c.active = true;
    const int base = lua_gettop(L);
    if (lua_checkstack(L, 4) == 0) {
        throw Error(Error::Kind::OutOfMemory, c.file, 0, std::string(kOutOfMemory));
    }
    c.executed = 0;
    c.budget = impl.options.limits.instructions;
    c.timeExceeded = false;
    c.outOfMemory = false;
    c.timeFile.clear();
    c.timeLine = 0;
    c.handled = false;
    ProtectedCall call;
    call.body = &body;
    lua_pushcfunction(L, messageHandler);
    lua_pushcfunction(L, trampoline);
    lua_pushlightuserdata(L, &call);
    const int status = lua_pcall(L, 1, 0, base + 1);
    std::string object;
    if (status != LUA_OK && lua_type(L, -1) == LUA_TSTRING) {
        object = lua_tostring(L, -1);
    }
    lua_settop(L, base);
    if (c.timeExceeded) {  // also when the error was swallowed (e.g. by a finalizer)
        const bool located = c.timeLine > 0;
        throw Error(Error::Kind::OutOfTime, located ? c.timeFile : c.file, c.timeLine, std::string(kOutOfTime));
    }
    if (status == LUA_ERRMEM || c.outOfMemory) {
        lua_pushcfunction(L, collectGarbage);  // give the memory back for the next call
        lua_pcall(L, 0, 0, 0);
        lua_settop(L, base);
        throw Error(Error::Kind::OutOfMemory, c.file, 0, std::string(kOutOfMemory));
    }
    if (status == LUA_ERRRUN && c.handled) {
        throw Error(Error::Kind::Runtime, c.errorFile, c.errorLine, c.errorMessage);
    }
    if (status == LUA_ERRERR) {
        throw Error(Error::Kind::Runtime, c.file, 0, "error in error handling: " + object);
    }
    if (status != LUA_OK) {  // a C++ exception crossed Lua's frames (a C function without guarded)
        throw Error(Error::Kind::Runtime, c.file, 0,
                    object.empty() ? std::string("internal error: an exception crossed the Lua interpreter") : object);
    }
    if (call.exception) {
        std::rethrow_exception(call.exception);
    }
}

std::vector<Ref> State::run(std::string_view source) {
    std::vector<Ref> results;
    protect([&](lua_State* L) {
        Control& c = m_impl->control;
        const int status = luaL_loadbufferx(L, source.data(), source.size(), m_impl->chunkName.c_str(), "t");
        if (status == LUA_ERRMEM) {
            throw Error(Error::Kind::OutOfMemory, c.file, 0, std::string(kOutOfMemory));
        }
        if (status != LUA_OK) {
            const char* text = lua_tostring(L, -1);
            std::string message = text != nullptr ? text : "?";
            int line = 0;
            std::string rest;
            if (splitLocation(message, c.chunkShortSrc, line, rest)) {
                message = std::move(rest);
            }
            throw Error(Error::Kind::Syntax, c.file, line, std::move(message));
        }
        lua_call(L, 0, LUA_MULTRET);
        const int count = lua_gettop(L);
        results.reserve(static_cast<std::size_t>(count));
        for (int i = 1; i <= count; ++i) {
            results.push_back(makeRef(L, i));
        }
    });
    return results;
}

std::vector<core::Value> State::call(const Ref& function, std::span<const core::Value> args) {
    std::vector<core::Value> results;
    protect([&](lua_State* L) {
        if (args.size() > 200 || lua_checkstack(L, static_cast<int>(args.size()) + 4) == 0) {
            throw std::invalid_argument("lua::State::call: too many arguments");
        }
        push(L, function);
        for (const core::Value& arg : args) {
            valueToLua(L, arg, Access::ReadOnly);
        }
        lua_call(L, static_cast<int>(args.size()), LUA_MULTRET);
        const int count = lua_gettop(L);
        results.reserve(static_cast<std::size_t>(count));
        for (int i = 1; i <= count; ++i) {
            results.push_back(luaToValue(L, i, count == 1 ? "the result" : "result " + std::to_string(i)));
        }
    });
    return results;
}

Ref State::field(const Ref& table, std::string_view key) {
    Ref result;
    protect([&](lua_State* L) {
        push(L, table);
        if (lua_type(L, -1) != LUA_TTABLE) {
            throw Error(Error::Kind::Runtime, file(), 0,
                        "attempt to index a " + std::string(luaL_typename(L, -1)) + " value (field '" + std::string(key) + "')");
        }
        lua_pushlstring(L, key.data(), key.size());
        lua_gettable(L, -2);
        result = makeRef(L, -1);
    });
    return result;
}

Ref State::global(std::string_view name) {
    Ref result;
    protect([&](lua_State* L) {
        lua_pushglobaltable(L);
        lua_pushlstring(L, name.data(), name.size());
        lua_rawget(L, -2);
        result = makeRef(L, -1);
    });
    return result;
}

core::Value State::value(const Ref& ref) {
    core::Value result;
    protect([&](lua_State* L) {
        push(L, ref);
        result = luaToValue(L, -1);
    });
    return result;
}

std::string State::typeName(const Ref& ref) {
    std::string result;
    protect([&](lua_State* L) {
        push(L, ref);
        result = luaL_typename(L, -1);
    });
    return result;
}

Ref State::makeRef(lua_State* L, int index) {
    luaL_checkstack(L, 3, nullptr);  // luaL_ref uses two more
    lua_pushvalue(L, index);
    const int ref = luaL_ref(L, LUA_REGISTRYINDEX);
    return Ref(std::static_pointer_cast<void>(m_impl), ref);
}

void State::push(lua_State* L, const Ref& ref) {
    const std::shared_ptr<void> owner = ref.m_owner.lock();
    if (owner != nullptr && owner.get() != static_cast<void*>(m_impl.get())) {
        throw std::invalid_argument("lua::State: a reference of another state");
    }
    luaL_checkstack(L, 1, nullptr);
    if (owner == nullptr || ref.m_ref < 0) {
        lua_pushnil(L);
    } else {
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref.m_ref);
    }
}

// --- Ref -----------------------------------------------------------------------------------------

Ref::Ref(Ref&& other) noexcept : m_owner(std::move(other.m_owner)), m_ref(std::exchange(other.m_ref, LUA_NOREF)) {}

Ref& Ref::operator=(Ref&& other) noexcept {
    if (this != &other) {
        reset();
        m_owner = std::move(other.m_owner);
        m_ref = std::exchange(other.m_ref, LUA_NOREF);
    }
    return *this;
}

Ref::~Ref() {
    reset();
}

Ref::operator bool() const noexcept {
    return m_ref >= 0 && !m_owner.expired();
}

void Ref::reset() noexcept {
    if (m_ref >= 0) {
        if (const std::shared_ptr<void> owner = m_owner.lock()) {
            luaL_unref(static_cast<State::Impl*>(owner.get())->L, LUA_REGISTRYINDEX, m_ref);
        }
    }
    m_owner.reset();
    m_ref = LUA_NOREF;
}

// --- raise ---------------------------------------------------------------------------------------

void raise(lua_State* L, std::string_view message) {
    luaL_where(L, 1);
    lua_pushlstring(L, message.data(), message.size());
    lua_concat(L, 2);
    lua_error(L);
    std::abort();  // not reached: lua_error throws
}

}  // namespace mdbridge::workflow::lua
