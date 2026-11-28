#pragma once
// The string-hash seed of a Lua state, chosen by the host instead of Lua's default (the clock and
// the addresses of the state, which differ from run to run).
//
// This header is included before third_party/lua/src/lstate.c (src/workflow/lua/CMakeLists.txt:
// -include / /FI) and by host.cpp. lstate.c computes the seed with `luai_makeseed(L)` unless that
// macro is defined, which is the hook Lua provides for it; the vendored file stays unmodified.
// lua_newstate stores the allocator's userdata (`G(L)->ud`) before it asks for the seed, so the host
// passes a SeedSource there (its allocator state derives from it).

namespace mdbridge::workflow::lua::detail {

struct SeedSource {
    unsigned int hashSeed = 0;
};

inline unsigned int hashSeedOf(const void* allocatorData) noexcept {
    return allocatorData != nullptr ? static_cast<const SeedSource*>(allocatorData)->hashSeed : 0u;
}

}  // namespace mdbridge::workflow::lua::detail

#define luai_makeseed(L) (::mdbridge::workflow::lua::detail::hashSeedOf(G(L)->ud))
