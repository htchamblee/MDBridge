# The MDATK commands of the command line that need no Qt (init, info, list, analyze, export-pdb,
# density info): linked into the `mdbridge` program and run in-process by the tests, so they build
# and are tested without Qt (MDBRIDGE_BUILD_APP=OFF). Also the interface's Qt-free part
# (project_view: the reference frame, the 3D view's scene, the figure preview).
FetchContent_MakeAvailable(CLI11)

add_library(mdbridge_commands STATIC
    commands.cpp
    commands.hpp
    overrides.cpp
    project_view.cpp
    project_view.hpp
    render_pipeline.cpp
    render_pipeline.hpp)
add_library(mdbridge::commands ALIAS mdbridge_commands)
target_include_directories(mdbridge_commands PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_commands PUBLIC mdbridge::core mdbridge::io mdbridge::analysis mdbridge::render_panels CLI11::CLI11)
target_compile_features(mdbridge_commands PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_commands)
