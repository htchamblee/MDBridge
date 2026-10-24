# core: data model, project files, step log, settings (no Qt).
FetchContent_MakeAvailable(tomlplusplus)

configure_file(version.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/include/mdbridge/version.hpp" @ONLY)

add_library(mdbridge_core STATIC
    clock.cpp
    hash.cpp
    json.cpp
    log.cpp
    naming.cpp
    paths.cpp
    project.cpp
    project_file.cpp
    project_io.cpp
    project_model.inc
    python.cpp
    settings.cpp
    step_log.cpp
    step_types.cpp
    step_types.inc
    text.cpp
    toml_writer.cpp
    unicode_tables.inc
    value.cpp)
add_library(mdbridge::core ALIAS mdbridge_core)
target_include_directories(mdbridge_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/include")
target_link_libraries(mdbridge_core PRIVATE tomlplusplus::tomlplusplus)
target_compile_features(mdbridge_core PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_core)
