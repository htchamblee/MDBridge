# The MDATK commands of the command line that need no Qt (demo, init, info, list, clean, analyze,
# export-pdb, density info / fit / simulate, pdb, build, dict; dict pubchem's transport comes from the
# program) and MDBridge's citations: linked into the `mdbridge` program and run in-process by the tests, so they build
# and are tested without Qt (MDBRIDGE_BUILD_APP=OFF). Also the interface's Qt-free part
# (project_view: the reference frame, the 3D view's scene, the figure preview).
FetchContent_MakeAvailable(CLI11)

add_library(mdbridge_commands STATIC
    chem_commands.cpp
    chem_commands.hpp
    citations.cpp
    citations.hpp
    command_support.cpp
    command_support.hpp
    commands.cpp
    commands.hpp
    demo.cpp
    demo.hpp
    density_pipeline.cpp
    density_pipeline.hpp
    overrides.cpp
    project_view.cpp
    project_view.hpp
    pubchem.cpp
    pubchem.hpp
    render_pipeline.cpp
    render_pipeline.hpp)
add_library(mdbridge::commands ALIAS mdbridge_commands)
target_include_directories(mdbridge_commands PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_commands PUBLIC mdbridge::core mdbridge::io mdbridge::analysis mdbridge::chem mdbridge::render_panels CLI11::CLI11)
target_compile_features(mdbridge_commands PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_commands)
# the demo system's coordinates rounded as numpy rounds them: no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties(demo.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(demo.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
