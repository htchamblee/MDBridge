#pragma once
// The Lua host of the workflow generators (docs/specs/m7-plan.md 4.4, 4.5): one sandboxed Lua 5.4
// state per generator file, with a memory cap and an instruction budget per call, errors that carry
// the file and line, and core::Value <-> Lua conversion (convert.hpp).
//
// Sandbox: the globals hold only `mdb`, `string`, `table`, `math`, `utf8` and the base functions
// without dofile, loadfile, print, warn, collectgarbage and rawset. `load` reads text only (mode "t",
// whatever the caller asks). `math.random` starts from a fixed seed and `math.randomseed` needs an
// explicit one (no seeding from the clock). `string.dump` is removed. The io, os, package, debug and
// coroutine libraries are not even compiled in. `pcall` and `xpcall` cannot catch running out of time
// or memory: those errors always reach the host (and xpcall's handler is not called for them).
// `setmetatable` refuses metatables with `__gc`: Lua runs finalizers with hooks off, where the
// instruction budget could not stop them.
//
// `mdb` holds `mdb.null` (core::Value's null inside arrays and tables) and `mdb.sorted_pairs(t)` (the
// keys of `t` in sorted order: numbers, then strings by bytes); the API of section 4.4 is added on top
// of this host (milestone 7 W1-B).
//
// A State is not thread-safe; different States may be used on different threads at the same time.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "mdbridge/core/value.hpp"

struct lua_State;

namespace mdbridge::workflow::lua {

struct Limits {
    // Bytes the state may hold at once (Lua's allocator; plan 4.5: 256 MiB).
    std::size_t memoryBytes = std::size_t{256} << 20;
    // Lua VM instructions per call from C++ (State::protect and the functions built on it),
    // counted in steps of kInstructionStep.
    std::uint64_t instructions = 100'000'000;
};

// The granularity of the instruction budget (the count hook fires every this many instructions).
inline constexpr int kInstructionStep = 1000;

// The string-hash seed of new states. Fixed, so that the order of `pairs` over string keys is the
// same in every run; tests run generators under two seeds to show that no output depends on it.
inline constexpr unsigned int kDefaultHashSeed = 0x4d444231u;  // "MDB1"

struct Options {
    Limits limits;
    unsigned int hashSeed = kDefaultHashSeed;
};

// The texts of the step's issue when a generator exceeds a limit (plan 4.5), after the plan's
// "internal error while checking this step: ".
inline constexpr std::string_view kOutOfTime = "the generator ran out of time";
inline constexpr std::string_view kOutOfMemory = "the generator ran out of memory";
// The prefix MDATK's plan puts before an exception escaping a step's check (plan.py:290).
inline constexpr std::string_view kIssuePrefix = "internal error while checking this step: ";

// An error of a generator file: a syntax error, an error raised while it ran, a limit reached, or a
// value it gave that has no core::Value.
class Error : public std::runtime_error {
public:
    enum class Kind { Syntax, Runtime, OutOfMemory, OutOfTime, Conversion };

    // what() is "file:line: message", or "file: message" when the line is not known (0).
    Error(Kind kind, std::string file, int line, std::string message);

    Kind kind() const noexcept { return m_kind; }
    // The generator file (the name the State was given) or, for a chunk loaded with `load`, that
    // chunk's name.
    const std::string& file() const noexcept { return m_file; }
    int line() const noexcept { return m_line; }
    // The message without the file and line.
    const std::string& message() const noexcept { return m_message; }

    // The issue the plan reports for the step: kIssuePrefix + kOutOfTime / kOutOfMemory for the
    // limits, kIssuePrefix + what() for every other error.
    std::string issueText() const;

private:
    Kind m_kind;
    std::string m_file;
    int m_line;
    std::string m_message;
};

class State;

// A Lua value kept for C++ (a registry reference). It may outlive its State; it is then empty.
class Ref {
public:
    Ref() = default;
    Ref(Ref&& other) noexcept;
    Ref& operator=(Ref&& other) noexcept;
    Ref(const Ref&) = delete;
    Ref& operator=(const Ref&) = delete;
    ~Ref();

    // True for a reference to a value other than nil.
    explicit operator bool() const noexcept;

private:
    friend class State;
    Ref(std::weak_ptr<void> owner, int ref) : m_owner(std::move(owner)), m_ref(ref) {}
    void reset() noexcept;

    std::weak_ptr<void> m_owner;
    int m_ref = -2;  // LUA_NOREF
};

class State {
public:
    // A new sandboxed state for the generator file `fileName` (a user step's path, or the name of a
    // built-in generator such as "steps/amber.lua"): the name its chunks and errors carry.
    explicit State(std::string fileName, Options options = {});
    ~State();
    State(const State&) = delete;
    State& operator=(const State&) = delete;

    const std::string& file() const noexcept;
    const Options& options() const noexcept;

    // Runs `source` (Lua text; binary chunks are refused) as the file's main chunk, with a fresh
    // instruction budget, and returns its results. Throws Error.
    std::vector<Ref> run(std::string_view source);

    // Calls `function` with `args` (read-only tables, convert.hpp) and converts its results.
    // Throws Error (also for results that have no core::Value).
    std::vector<core::Value> call(const Ref& function, std::span<const core::Value> args);

    // The field `key` of the table (or read-only table) `table`; empty if absent. Throws Error.
    Ref field(const Ref& table, std::string_view key);
    // The global `name` of the state; empty if absent.
    Ref global(std::string_view name);
    // `ref` converted to a core::Value. Throws Error.
    core::Value value(const Ref& ref);
    // The Lua type name of `ref`'s value ("nil", "table", "function", ...).
    std::string typeName(const Ref& ref);

    // Runs `body` in protected mode with a fresh instruction budget: Lua errors raised in `body`
    // (including in Lua functions it calls with lua_call) and the limits become Error, C++
    // exceptions thrown by `body` propagate unchanged. The stack is restored afterwards. This is the
    // way to use the raw state: the functions below and in convert.hpp need it.
    void protect(const std::function<void(lua_State*)>& body);

    // Inside `protect`: a reference to the value at `index`, and pushing a reference's value.
    Ref makeRef(lua_State* L, int index);
    void push(lua_State* L, const Ref& ref);

    // Bytes the state holds now.
    std::size_t memoryInUse() const noexcept;

    struct Impl;

private:
    std::shared_ptr<Impl> m_impl;
};

// Raises a Lua error with `message`, prefixed with the file and line of the calling Lua function
// (as luaL_error does). Never returns.
[[noreturn]] void raise(lua_State* L, std::string_view message);

// The body of a lua_CFunction written in C++: C++ exceptions (std::exception) thrown by `body` become
// Lua errors with their what() text, raised after the exception is destroyed. Lua's own errors pass
// through (they are not std::exception). Use as `return guarded(L, [&] { ...; return results; });`.
template <class Body>
int guarded(lua_State* L, Body&& body) {
    std::string message;
    try {
        return std::forward<Body>(body)();
    } catch (const std::exception& error) {
        message = error.what();
    }
    raise(L, message);
}

}  // namespace mdbridge::workflow::lua
