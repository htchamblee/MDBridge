# workflow/lua: the Lua host of milestone 7's step generators (docs/specs/m7-plan.md 4.4, 4.5): the
# vendored Lua 5.4 interpreter (mdbridge_lua) and MDBridge's host library on it (mdbridge_lua_host:
# the sandboxed state, its memory and instruction limits, errors with file and line, core::Value <->
# Lua conversion, the determinism lint). No Qt.

# Lua 5.4.9 as PUC-Rio releases it (tools/vendor_lua.py): third-party code, never edited; every file is
# checked against SHA256SUMS here, so an edited copy stops the configuration
set(MDBRIDGE_LUA_DIR "${PROJECT_SOURCE_DIR}/third_party/lua")
file(STRINGS "${MDBRIDGE_LUA_DIR}/SHA256SUMS" lua_sums)
foreach(line IN LISTS lua_sums)
    string(REGEX MATCH "^([0-9a-f]+)  (.+)$" matched "${line}")
    file(SHA256 "${MDBRIDGE_LUA_DIR}/${CMAKE_MATCH_2}" actual)
    if(NOT actual STREQUAL CMAKE_MATCH_1)
        message(FATAL_ERROR "third_party/lua/${CMAKE_MATCH_2} differs from Lua 5.4.9's copy "
                            "(SHA256SUMS; recreate it with tools/vendor_lua.py)")
    endif()
endforeach()

# The interpreter and the libraries the sandbox opens (base, string, table, math, utf8), compiled as
# C++ so that Lua errors unwind with exceptions and the destructors of the host's C++ frames run.
# Left out: the io, os, package (loadlib), debug and coroutine libraries and linit.c, so no generator
# can reach files, programs, the clock or the environment even through a bug in the host.
set(lua_sources
    lapi.c lcode.c lctype.c ldebug.c ldo.c ldump.c lfunc.c lgc.c llex.c lmem.c lobject.c lopcodes.c
    lparser.c lstate.c lstring.c ltable.c ltm.c lundump.c lvm.c lzio.c
    lauxlib.c lbaselib.c lmathlib.c lstrlib.c ltablib.c lutf8lib.c)
list(TRANSFORM lua_sources PREPEND "${MDBRIDGE_LUA_DIR}/src/")
set_source_files_properties(${lua_sources} PROPERTIES LANGUAGE CXX)  # CMP0119: -x c++ / /TP
add_library(mdbridge_lua STATIC ${lua_sources})
# SYSTEM: MDBridge's code includes lua.h with <...>, so its strict warnings (and MSVC's /external:W0)
# do not apply to Lua's headers
target_include_directories(mdbridge_lua SYSTEM PUBLIC "${MDBRIDGE_LUA_DIR}/src")
target_compile_features(mdbridge_lua PUBLIC cxx_std_20)
set_target_properties(mdbridge_lua PROPERTIES POSITION_INDEPENDENT_CODE ON)
# The string-hash seed comes from the host (Options::hashSeed) instead of the clock and addresses:
# lstate.c takes luai_makeseed from lua_seed.hpp, which is included before it (Lua's own hook,
# "#if !defined(luai_makeseed)"; the file stays unmodified)
if(MSVC)
    set_source_files_properties("${MDBRIDGE_LUA_DIR}/src/lstate.c" PROPERTIES
        COMPILE_OPTIONS "/FI${CMAKE_CURRENT_SOURCE_DIR}/lua_seed.hpp")
else()
    set_source_files_properties("${MDBRIDGE_LUA_DIR}/src/lstate.c" PROPERTIES
        COMPILE_OPTIONS "-include;${CMAKE_CURRENT_SOURCE_DIR}/lua_seed.hpp")
endif()
# Warnings: GCC 11 compiles these sources as C++ without a warning at -Wall -Wextra, which stay on
# (not -Werror: the code is not MDBridge's to fix). MDBridge's C++ style warnings (-Wold-style-cast,
# -Wpedantic for the computed gotos of lvm.c, ...) are not applied: the sources are C. MSVC keeps its
# default level /W3, without /W4 and /WX, for the same reason.
if(NOT MSVC)
    target_compile_options(mdbridge_lua PRIVATE -Wall -Wextra)
endif()
# Debug builds check the use of Lua's C API (lapi.c's api_check) by the host
target_compile_definitions(mdbridge_lua PRIVATE $<$<CONFIG:Debug>:LUA_USE_APICHECK>)

add_library(mdbridge_lua_host STATIC
    convert.cpp
    host.cpp
    host_internal.hpp
    include/mdbridge/workflow/lua/convert.hpp
    include/mdbridge/workflow/lua/host.hpp
    include/mdbridge/workflow/lua/lint.hpp
    lint.cpp
    lua_seed.hpp)
add_library(mdbridge::lua_host ALIAS mdbridge_lua_host)
target_include_directories(mdbridge_lua_host PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_lua_host PUBLIC mdbridge::core mdbridge_lua)
target_compile_features(mdbridge_lua_host PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_lua_host)

# The determinism lint as a program: `mdbridge_lua_lint FILE_OR_FOLDER...` prints each finding as
# "file:line:column: rule: message" and exits with 1 if there is any (ctest runs it over the bundled
# generators, src/workflow/lua)
if(MDBRIDGE_BUILD_TESTS)
    add_executable(mdbridge_lua_lint lint_main.cpp)
    target_link_libraries(mdbridge_lua_lint PRIVATE mdbridge::lua_host)
    mdbridge_target_warnings(mdbridge_lua_lint)
endif()
