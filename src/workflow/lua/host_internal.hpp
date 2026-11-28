#pragma once
// Shared by host.cpp and convert.cpp.

#include <string>

struct lua_State;

namespace mdbridge::workflow::lua::detail {

// The file name the state was created for (State::file()).
const std::string& stateFile(lua_State* L);

}  // namespace mdbridge::workflow::lua::detail
