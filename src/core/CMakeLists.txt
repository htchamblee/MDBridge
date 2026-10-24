# core: data model, project files, step log, settings (no Qt).
FetchContent_MakeAvailable(tomlplusplus)

configure_file(version.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/include/mdbridge/version.hpp" @ONLY)

add_library(mdbridge_core STATIC
    json.cpp
    paths.cpp
    project_file.cpp
    value.cpp)
add_library(mdbridge::core ALIAS mdbridge_core)
target_include_directories(mdbridge_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/include")
target_link_libraries(mdbridge_core PRIVATE tomlplusplus::tomlplusplus)
target_compile_features(mdbridge_core PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_core)
