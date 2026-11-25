set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)  # MSVC: the same runtime as MDBridge
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# core: project files read exactly as MDATK reads them (fixtures in data/, see data/README.md)
add_executable(mdbridge_core_tests
    core/test_atom_selection.cpp
    core/test_citations.cpp
    core/test_conda.cpp
    core/test_dataset.cpp
    core/test_hash.cpp
    core/test_json.cpp
    core/test_naming.cpp
    core/test_ntpath.cpp
    core/test_panels.cpp
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
    io/test_cif.cpp
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

# chem: secondary structure (dot-bracket, builtin base pairs) as MDATK computes it
add_executable(mdbridge_chem_tests chem/test_secstruct.cpp)
target_include_directories(mdbridge_chem_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_chem_tests PRIVATE mdbridge::chem GTest::gtest_main)
target_compile_definitions(mdbridge_chem_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_chem_tests)
gtest_discover_tests(mdbridge_chem_tests DISCOVERY_MODE PRE_TEST)


# analysis: the framework (registry, instance ids, labels, cache keys, runner, common helpers) as
# MDATK computes it (data/analysis, tests/parity/make_analysis_framework_cases.py), the analyses
# against MDATK's runs (tests/parity/make_analysis_cases.py), MDATK's spatial helpers
# (data/analysis/spatial_extra_cases.json), the periodic pair search against MDAnalysis
# (data/analysis/distance_cases.json) and the optimal-superposition RMSD against MDAnalysis
# (data/analysis/rmsd_cases.json, docs/specs/m3-rmsd.md)
add_executable(mdbridge_analysis_tests
    analysis/test_basepairs.cpp
    analysis/test_contacts_hbonds.cpp
    analysis/test_energy.cpp
    analysis/test_framework.cpp
    analysis/test_improvements_analyses.cpp
    analysis/test_mg_stacking.cpp
    analysis/test_parallel.cpp
    analysis/test_periodic.cpp
    analysis/test_pipeline.cpp
    analysis/test_rmsd_analysis.cpp
    analysis/test_rmsd_qcp.cpp
    analysis/test_spatial.cpp
    analysis/test_structure.cpp)
target_include_directories(mdbridge_analysis_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_analysis_tests PRIVATE mdbridge::analysis mdbridge::io GTest::gtest_main)
target_compile_definitions(mdbridge_analysis_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_analysis_tests)
# the spatial tests compute their expected sums and products as numpy does: no fused multiply-adds
if(MSVC)
    set_source_files_properties(analysis/test_spatial.cpp analysis/test_periodic.cpp
                                PROPERTIES COMPILE_OPTIONS "/fp:precise")
else()
    set_source_files_properties(analysis/test_spatial.cpp analysis/test_periodic.cpp
                                PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
endif()
# every case twice: with the frame-parallel runner (4 threads, whatever the machine) and
# sequentially (MDBRIDGE_THREADS=1, tests named "<test>/threads1"); the numerics that do not go
# through the runner (spatial helpers, periodic pair search, QCP RMSD) only once
gtest_discover_tests(mdbridge_analysis_tests DISCOVERY_MODE PRE_TEST PROPERTIES ENVIRONMENT "MDBRIDGE_THREADS=4")
set(mdbridge_runner_free_suites
    Kabsch Linalg MinimumImageVectors NumpySum PairsWithin PairsWithinPeriodic UsableBox RmsdQcp)
list(TRANSFORM mdbridge_runner_free_suites APPEND ".*")
list(JOIN mdbridge_runner_free_suites ":" mdbridge_runner_free_filter)
gtest_discover_tests(mdbridge_analysis_tests DISCOVERY_MODE PRE_TEST TEST_SUFFIX "/threads1"
                     TEST_FILTER "-${mdbridge_runner_free_filter}"
                     PROPERTIES ENVIRONMENT "MDBRIDGE_THREADS=1")

# app/commands: MDATK's init, info, list, analyze, export-pdb and density info, run in-process
# against what MDATK printed and wrote (data/cli, tests/parity/cli_parity.py); no Qt needed
add_executable(mdbridge_cli_tests cli/test_commands.cpp)
target_include_directories(mdbridge_cli_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_cli_tests PRIVATE mdbridge::commands GTest::gtest_main)
target_compile_definitions(mdbridge_cli_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_cli_tests)
gtest_discover_tests(mdbridge_cli_tests DISCOVERY_MODE PRE_TEST)

# the same commands as a program, for tests/parity/cli_parity.py --mdbridge on builds without Qt
add_executable(mdbridge_commands_cli cli/commands_cli.cpp)
target_link_libraries(mdbridge_commands_cli PRIVATE mdbridge::commands)
mdbridge_target_warnings(mdbridge_commands_cli)

# render core (no Qt): colours, colormaps, chain palettes and themes (data/render/colors,
# tests/parity/make_color_cases.py); axis limits, ticks and tick labels (data/render/ticks,
# tests/parity/make_tick_cases.py); the bundled fonts and text metrics (data/render/text,
# tests/parity/make_text_cases.py; the fonts in data/fonts), against Matplotlib 3.11.2 and MDATK
add_executable(mdbridge_render_core_tests render/test_colors.cpp render/test_text.cpp render/test_ticks.cpp)
target_include_directories(mdbridge_render_core_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_render_core_tests PRIVATE mdbridge::render_core mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_render_core_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_render_core_tests)
gtest_discover_tests(mdbridge_render_core_tests DISCOVERY_MODE PRE_TEST)

# render/molecule: the 3D representations' geometry (topology roles, styles, Gaussian surface,
# scene colours, the scipy filters) against MDATK (data/render/molecule_geometry,
# tests/parity/make_molecule_geometry_cases.py); no Qt
add_executable(mdbridge_molecule_tests render/test_molecule_geometry.cpp)
target_link_libraries(mdbridge_molecule_tests PRIVATE mdbridge::render_molecule mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_molecule_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_molecule_tests)
gtest_discover_tests(mdbridge_molecule_tests DISCOVERY_MODE PRE_TEST)

# render/video: MP4 encoding with the LGPL FFmpeg libraries (codec and quality map, LGPL-only check,
# encoded files checked by demuxing); the encoding tests skip in builds without FFmpeg. No Qt.
add_executable(mdbridge_video_tests render/test_video_encoder.cpp)
target_link_libraries(mdbridge_video_tests PRIVATE mdbridge::video GTest::gtest_main)
mdbridge_target_warnings(mdbridge_video_tests)
if(WIN32 AND MDBRIDGE_VIDEO_RUNTIME_FILES)  # the FFmpeg DLLs next to the test program
    add_custom_command(TARGET mdbridge_video_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MDBRIDGE_VIDEO_RUNTIME_FILES} $<TARGET_FILE_DIR:mdbridge_video_tests>
        VERBATIM)
endif()
gtest_discover_tests(mdbridge_video_tests DISCOVERY_MODE PRE_TEST)

# the committed SPIR-V matches the shader sources (tools/compile_shaders.py)
add_test(NAME render.shaders_up_to_date
         COMMAND ${CMAKE_COMMAND} -DSHADER_DIR=${PROJECT_SOURCE_DIR}/src/render/shaders
                 -P ${PROJECT_SOURCE_DIR}/cmake/CheckShaders.cmake)

if(MDBRIDGE_BUILD_APP)
    add_executable(mdbridge_render_tests render/test_vulkan_probe.cpp)
    target_link_libraries(mdbridge_render_tests PRIVATE mdbridge::render Vulkan::Headers GTest::gtest_main)
    mdbridge_target_warnings(mdbridge_render_tests)
    gtest_discover_tests(mdbridge_render_tests DISCOVERY_MODE PRE_TEST)

    # the figure fonts in Qt: resources, registration, text drawn from outlines, the batch
    # application on the offscreen platform (no display needed)
    add_executable(mdbridge_text_qt_tests render/test_text_qt.cpp)
    target_link_libraries(mdbridge_text_qt_tests PRIVATE mdbridge::render GTest::gtest)
    mdbridge_target_warnings(mdbridge_text_qt_tests)
    gtest_discover_tests(mdbridge_text_qt_tests DISCOVERY_MODE PRE_TEST)
endif()
