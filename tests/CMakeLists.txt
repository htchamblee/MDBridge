set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)  # MSVC: the same runtime as MDBridge
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# core: project files read exactly as MDATK reads them (fixtures in data/, see data/README.md)
add_executable(mdbridge_core_tests
    core/test_atom_selection.cpp
    core/test_dataset.cpp
    core/test_hash.cpp
    core/test_json.cpp
    core/test_naming.cpp
    core/test_project_file.cpp
    core/test_project_io.cpp
    core/test_progress.cpp
    core/test_project_model.cpp
    core/test_python.cpp
    core/test_rigid_transform.cpp
    core/test_selections.cpp
    core/test_step_log.cpp
    core/test_step_types.cpp
    core/test_text.cpp
    core/test_toml_writer.cpp
    core/test_topology.cpp)
target_link_libraries(mdbridge_core_tests PRIVATE mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_core_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_core_tests)
gtest_discover_tests(mdbridge_core_tests DISCOVERY_MODE PRE_TEST)

# io: topologies, trajectories, PDB files and tables read and written as MDATK does
add_executable(mdbridge_io_tests
    io/test_external.cpp
    io/test_loader.cpp
    io/test_mrc.cpp
    io/test_netcdf.cpp
    io/test_parm7.cpp
    io/test_pdb.cpp
    io/test_pdb_export.cpp
    io/test_trajectory.cpp)
target_include_directories(mdbridge_io_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_io_tests PRIVATE mdbridge::io GTest::gtest_main)
target_compile_definitions(mdbridge_io_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_io_tests)
gtest_discover_tests(mdbridge_io_tests DISCOVERY_MODE PRE_TEST)

# the committed SPIR-V matches the shader sources (tools/compile_shaders.py)
add_test(NAME render.shaders_up_to_date
         COMMAND ${CMAKE_COMMAND} -DSHADER_DIR=${PROJECT_SOURCE_DIR}/src/render/shaders
                 -P ${PROJECT_SOURCE_DIR}/cmake/CheckShaders.cmake)

if(MDBRIDGE_BUILD_APP)
    add_executable(mdbridge_render_tests render/test_vulkan_probe.cpp)
    target_link_libraries(mdbridge_render_tests PRIVATE mdbridge::render Vulkan::Headers GTest::gtest_main)
    mdbridge_target_warnings(mdbridge_render_tests)
    gtest_discover_tests(mdbridge_render_tests DISCOVERY_MODE PRE_TEST)
endif()
