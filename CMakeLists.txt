# MDBridge: C++ successor of MDATK (see docs/ARCHITECTURE.md).
#
#   cmake -S . -B build -G Ninja -DCMAKE_PREFIX_PATH=<Qt 6>/lib/cmake -DCMAKE_BUILD_TYPE=Release
#   cmake --build build && ctest --test-dir build
#
# Qt 6 must be installed; every other dependency is fetched at configure time (cmake/Dependencies.cmake).
# docs/BUILDING.md has the details.
cmake_minimum_required(VERSION 3.24)

project(MDBridge
    VERSION 0.1.0
    DESCRIPTION "Analysis and visualisation of MD trajectories of RNA systems (successor of MDATK)"
    HOMEPAGE_URL "https://github.com/htchamblee/MDBridge"
    LANGUAGES CXX)

# The MDATK release whose files and results MDBridge must reproduce.
set(MDBRIDGE_MDATK_REFERENCE "0.0.4.1")

option(MDBRIDGE_BUILD_APP "Build the Qt front ends (GUI and CLI); OFF builds the core libraries only" ON)
option(MDBRIDGE_BUILD_TESTS "Build the unit tests" ON)
option(MDBRIDGE_USE_SYSTEM_DEPS "Prefer installed toml++, CLI11, GoogleTest and Vulkan headers over fetching them" OFF)
option(MDBRIDGE_WARNINGS_AS_ERRORS "Treat compiler warnings in MDBridge's own code as errors" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(CompilerOptions)
include(Dependencies)

# Linux installs are relocatable (AppImage): programs find the bundled libraries next to them.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")

if(MDBRIDGE_BUILD_TESTS)
    enable_testing()
endif()

add_subdirectory(src/core)
if(MDBRIDGE_BUILD_APP)
    find_package(Qt6 6.5 REQUIRED COMPONENTS Core Gui Widgets)
    add_subdirectory(src/render)
    add_subdirectory(src/app)
endif()
if(MDBRIDGE_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# licence texts: MDBridge's notice, Qt's LGPL, and those of the fetched dependencies compiled in
if(WIN32)  # next to bin\, as .txt files that open in Notepad
    set(MDBRIDGE_DOC_DIR .)
    install(FILES LICENSE DESTINATION . RENAME LICENSE.txt)
    install(FILES THIRD_PARTY_NOTICES.md DESTINATION . RENAME THIRD_PARTY_NOTICES.txt)
else()
    set(MDBRIDGE_DOC_DIR ${CMAKE_INSTALL_DATADIR}/doc/mdbridge)
    install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION ${MDBRIDGE_DOC_DIR})
endif()
install(DIRECTORY packaging/licenses/ DESTINATION ${MDBRIDGE_DOC_DIR}/licenses)
foreach(dep tomlplusplus CLI11 VulkanHeaders)
    FetchContent_GetProperties(${dep} POPULATED populated SOURCE_DIR source)
    if(populated)  # not when installed headers or packages were used instead
        file(GLOB texts "${source}/LICENSE*" "${source}/LICENSES/*")
        foreach(text IN LISTS texts)
            if(NOT IS_DIRECTORY "${text}")
                get_filename_component(name "${text}" NAME_WLE)
                string(REPLACE "VulkanHeaders" "Vulkan-Headers" label "${dep}")
                install(FILES "${text}" DESTINATION ${MDBRIDGE_DOC_DIR}/licenses RENAME "${label}-${name}.txt")
            endif()
        endforeach()
    endif()
endforeach()
