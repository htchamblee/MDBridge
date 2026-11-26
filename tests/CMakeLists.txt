set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)  # MSVC: the same runtime as MDBridge
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
if(MDBRIDGE_BUILD_APP)
    # Qt Test, once, for every GUI test executable (mdbridge_gui_*_tests: label gui, offscreen platform)
    find_package(Qt6 REQUIRED COMPONENTS Test)
endif()

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
    core/test_project_document.cpp
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

# app/commands: MDATK's demo, init, info, list, clean, analyze, export-pdb and density info, run
# in-process against what MDATK printed and wrote (data/cli, tests/parity/cli_parity.py), and the demo,
# cleaning, process runner and citations on their own (cli/test_cleaning.cpp); no Qt needed
add_executable(mdbridge_cli_tests cli/test_commands.cpp cli/test_cleaning.cpp)
target_include_directories(mdbridge_cli_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_cli_tests PRIVATE mdbridge::commands GTest::gtest_main)
target_compile_definitions(mdbridge_cli_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_cli_tests)
gtest_discover_tests(mdbridge_cli_tests DISCOVERY_MODE PRE_TEST)

# render, charts and preview with a recording backend (the drawing is the Qt program's: cli_parity.py)
# and the interface's project data without Qt (src/app/commands/project_view.cpp)
add_executable(mdbridge_render_cli_tests cli/test_render_commands.cpp cli/test_project_view.cpp cli/test_demo_commands.cpp)
target_include_directories(mdbridge_render_cli_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_render_cli_tests PRIVATE mdbridge::commands GTest::gtest_main)
target_compile_definitions(mdbridge_render_cli_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_render_cli_tests)
gtest_discover_tests(mdbridge_render_cli_tests DISCOVERY_MODE PRE_TEST)

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

# render core: the figure framework (model, constrained layout, record writer, display list) on the
# chrome of MDATK's figures (data/render/figures, tests/parity/make_figure_cases.py); no Qt
add_executable(mdbridge_figure_tests render/test_figure.cpp)
target_include_directories(mdbridge_figure_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_figure_tests PRIVATE mdbridge::render_core mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_figure_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_figure_tests)
gtest_discover_tests(mdbridge_figure_tests DISCOVERY_MODE PRE_TEST)

# render/panels: the figure's panels and build_figure against MDATK's figures (data/render/figures,
# tests/parity/make_figure_cases.py; figure_records.hpp reads the cases and their inputs); no Qt
add_executable(mdbridge_panel_tests render/test_plot_panels.cpp)
target_include_directories(mdbridge_panel_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_panel_tests PRIVATE mdbridge::render_panels mdbridge::analysis mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_panel_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_panel_tests)
gtest_discover_tests(mdbridge_panel_tests DISCOVERY_MODE PRE_TEST)

# render/panels: the sequence and secondary-structure panels (counterparts of MDATK's
# tests/test_seqpanels.py; contents against data/render/figures, tests/parity/make_figure_cases.py); no Qt
add_executable(mdbridge_seqpanel_tests render/test_seqpanels.cpp)
target_include_directories(mdbridge_seqpanel_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_seqpanel_tests PRIVATE mdbridge::render_panels mdbridge::analysis mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_seqpanel_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_seqpanel_tests)
gtest_discover_tests(mdbridge_seqpanel_tests DISCOVERY_MODE PRE_TEST)

# render/panels: the standalone charts against MDATK's render_charts (data/render/charts,
# tests/parity/make_chart_cases.py): paths, files, log, progress, warnings, figures; no Qt
add_executable(mdbridge_chart_tests render/test_charts.cpp)
target_include_directories(mdbridge_chart_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_chart_tests PRIVATE mdbridge::render_panels mdbridge::analysis mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_chart_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_chart_tests)
gtest_discover_tests(mdbridge_chart_tests DISCOVERY_MODE PRE_TEST)

# render core: residue shading, orientation, the example dataset and numpy's PCG64 against MDATK and
# numpy (data/render/shading, tests/parity/make_shading_cases.py); no Qt
add_executable(mdbridge_shading_tests render/test_shading.cpp)
target_link_libraries(mdbridge_shading_tests PRIVATE mdbridge::render_core mdbridge::analysis mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_shading_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_shading_tests)
gtest_discover_tests(mdbridge_shading_tests DISCOVERY_MODE PRE_TEST)

# render/molecule: the 3D representations' geometry (topology roles, styles, Gaussian surface,
# scene colours, the scipy filters) against MDATK (data/render/molecule_geometry,
# tests/parity/make_molecule_geometry_cases.py); no Qt
add_executable(mdbridge_molecule_tests render/test_molecule_geometry.cpp)
target_link_libraries(mdbridge_molecule_tests PRIVATE mdbridge::render_molecule mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_molecule_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_molecule_tests)
gtest_discover_tests(mdbridge_molecule_tests DISCOVERY_MODE PRE_TEST)

# render/molecule: the software 3D renderer (counterparts of MDATK's tests/test_molecule_render.py;
# RGBA against MDATK's renders, data/render/molecule, tests/parity/make_molecule_cases.py); no Qt
add_executable(mdbridge_molecule_render_tests render/test_molecule_render.cpp)
target_link_libraries(mdbridge_molecule_render_tests PRIVATE mdbridge::render_molecule mdbridge::render_core
                      mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_molecule_render_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_molecule_render_tests)
gtest_discover_tests(mdbridge_molecule_render_tests DISCOVERY_MODE PRE_TEST)

# render/panels: the structure view on the figure framework against MDATK's figures
# (data/render/figures, tests/parity/make_figure_cases.py): 2D trace, tags, legend, density shadow,
# 3D layer (software renderer); no Qt
add_executable(mdbridge_structure_panel_tests render/test_structure_panel.cpp)
target_include_directories(mdbridge_structure_panel_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_structure_panel_tests PRIVATE mdbridge::render_panels mdbridge::analysis mdbridge::core
                      GTest::gtest_main)
target_compile_definitions(mdbridge_structure_panel_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_structure_panel_tests)
gtest_discover_tests(mdbridge_structure_panel_tests DISCOVERY_MODE PRE_TEST)

# render/molecule: the density shadow (project_density, shadow_rgba, anchored, scipy's
# affine_transform) against MDATK (data/render/density_shadow,
# tests/parity/make_density_shadow_cases.py); no Qt
add_executable(mdbridge_density_shadow_tests render/test_density_shadow.cpp)
target_link_libraries(mdbridge_density_shadow_tests PRIVATE mdbridge::render_molecule mdbridge::io mdbridge::core
                      GTest::gtest_main)
target_compile_definitions(mdbridge_density_shadow_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_density_shadow_tests)
gtest_discover_tests(mdbridge_density_shadow_tests DISCOVERY_MODE PRE_TEST)

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

# render/video: the frame loop (MDATK's render_video) with a test compositor: positions, order,
# threads, progress, cancellation, whole runs of data/render/cli's video cases through the encoder
# (tests/parity/make_render_cli_cases.py). No Qt.
add_executable(mdbridge_video_loop_tests render/test_video_loop.cpp)
target_include_directories(mdbridge_video_loop_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
target_link_libraries(mdbridge_video_loop_tests PRIVATE mdbridge::video_loop mdbridge::commands mdbridge::analysis
                      mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_video_loop_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                           MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
mdbridge_target_warnings(mdbridge_video_loop_tests)
if(WIN32 AND MDBRIDGE_VIDEO_RUNTIME_FILES)  # the FFmpeg DLLs next to the test program
    add_custom_command(TARGET mdbridge_video_loop_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MDBRIDGE_VIDEO_RUNTIME_FILES} $<TARGET_FILE_DIR:mdbridge_video_loop_tests>
        VERBATIM)
endif()
gtest_discover_tests(mdbridge_video_loop_tests DISCOVERY_MODE PRE_TEST)

# render/vulkan: Vulkan without a window (loader, device, offscreen target and read-back, pipelines,
# pipeline cache) and the probe; no Qt. The GPU tests run on the machine's Vulkan 1.3 driver
# (lavapipe in CI) and skip without one, unless MDBRIDGE_REQUIRE_VULKAN=1 (tests/render/vulkan_harness.hpp)
add_executable(mdbridge_vulkan_tests render/test_vulkan_offscreen.cpp)
target_include_directories(mdbridge_vulkan_tests PRIVATE "${PROJECT_SOURCE_DIR}/src/render")  # shaders/spirv.hpp
target_link_libraries(mdbridge_vulkan_tests PRIVATE mdbridge::render_vulkan GTest::gtest_main)
mdbridge_target_warnings(mdbridge_vulkan_tests)
gtest_discover_tests(mdbridge_vulkan_tests DISCOVERY_MODE PRE_TEST PROPERTIES LABELS vulkan)

# render/vulkan: the molecule scene on the GPU against the software renderer (data/render/molecule:
# coverage IoU >= 0.99, SSIM >= 0.97), the GPU stage renderer and its renderer choice; no Qt. On
# lavapipe in CI (label vulkan); MDBRIDGE_VULKAN_BENCHMARK=1 prints frames/s
add_executable(mdbridge_vulkan_scene_tests render/test_vulkan_scene.cpp)
target_link_libraries(mdbridge_vulkan_scene_tests PRIVATE mdbridge::render_panels_vulkan mdbridge::render_core
                      mdbridge::core GTest::gtest_main)
target_compile_definitions(mdbridge_vulkan_scene_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
mdbridge_target_warnings(mdbridge_vulkan_scene_tests)
gtest_discover_tests(mdbridge_vulkan_scene_tests DISCOVERY_MODE PRE_TEST PROPERTIES LABELS vulkan)

# the committed SPIR-V matches the shader sources (tools/compile_shaders.py)
add_test(NAME render.shaders_up_to_date
         COMMAND ${CMAKE_COMMAND} -DSHADER_DIR=${PROJECT_SOURCE_DIR}/src/render/shaders
                 -P ${PROJECT_SOURCE_DIR}/cmake/CheckShaders.cmake)

if(MDBRIDGE_BUILD_APP)
    add_executable(mdbridge_render_tests render/test_vulkan_probe.cpp)
    target_link_libraries(mdbridge_render_tests PRIVATE mdbridge::render Vulkan::Headers GTest::gtest_main)
    mdbridge_target_warnings(mdbridge_render_tests)
    gtest_discover_tests(mdbridge_render_tests DISCOVERY_MODE PRE_TEST)

    # the GUI's Vulkan 3D view drawing a molecule (GPU scene renderer on Qt's device); needs a
    # display with Vulkan (xvfb-run + lavapipe), skips without one
    add_executable(mdbridge_vulkan_window_tests render/test_vulkan_window.cpp)
    target_link_libraries(mdbridge_vulkan_window_tests PRIVATE mdbridge::render GTest::gtest)
    mdbridge_target_warnings(mdbridge_vulkan_window_tests)
    gtest_discover_tests(mdbridge_vulkan_window_tests DISCOVERY_MODE PRE_TEST)

    # the figure fonts in Qt: resources, registration, text drawn from outlines, the batch
    # application on the offscreen platform (no display needed)
    add_executable(mdbridge_text_qt_tests render/test_text_qt.cpp)
    target_link_libraries(mdbridge_text_qt_tests PRIVATE mdbridge::render GTest::gtest)
    mdbridge_target_warnings(mdbridge_text_qt_tests)
    gtest_discover_tests(mdbridge_text_qt_tests DISCOVERY_MODE PRE_TEST)

    # the figure framework's QPainter backends: raster images of exact size, SVG with text as text,
    # PDF pages of the figure's size, the frame loop's layers (offscreen platform)
    add_executable(mdbridge_figure_qt_tests render/test_figure_qt.cpp)
    target_link_libraries(mdbridge_figure_qt_tests PRIVATE mdbridge::render GTest::gtest)
    target_compile_definitions(mdbridge_figure_qt_tests PRIVATE MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
    mdbridge_target_warnings(mdbridge_figure_qt_tests)
    gtest_discover_tests(mdbridge_figure_qt_tests DISCOVERY_MODE PRE_TEST)

    # the standalone charts' files from the Qt backend: every format MDATK's charts write, their
    # first bytes and sizes as MDATK's (data/render/charts), PostScript and TIFF (offscreen platform)
    add_executable(mdbridge_charts_qt_tests render/test_charts_qt.cpp)
    target_include_directories(mdbridge_charts_qt_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
    target_link_libraries(mdbridge_charts_qt_tests PRIVATE mdbridge::render mdbridge::render_panels mdbridge::analysis GTest::gtest)
    target_compile_definitions(mdbridge_charts_qt_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                               MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
    mdbridge_target_warnings(mdbridge_charts_qt_tests)
    gtest_discover_tests(mdbridge_charts_qt_tests DISCOVERY_MODE PRE_TEST)

    # the video frame loop with the QPainter compositor: MDATK's full-run video (24 frames at
    # 320 x 180), frames identical for any thread count, decoded frames vs composed (PSNR)
    add_executable(mdbridge_video_qt_tests render/test_video_qt.cpp)
    target_include_directories(mdbridge_video_qt_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
    target_link_libraries(mdbridge_video_qt_tests PRIVATE mdbridge::render mdbridge::commands mdbridge::analysis GTest::gtest)
    target_compile_definitions(mdbridge_video_qt_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                               MDBRIDGE_FONT_DIR="${PROJECT_SOURCE_DIR}/data/fonts")
    mdbridge_target_warnings(mdbridge_video_qt_tests)
    if(WIN32 AND MDBRIDGE_VIDEO_RUNTIME_FILES)
        add_custom_command(TARGET mdbridge_video_qt_tests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MDBRIDGE_VIDEO_RUNTIME_FILES} $<TARGET_FILE_DIR:mdbridge_video_qt_tests>
            VERBATIM)
    endif()
    gtest_discover_tests(mdbridge_video_qt_tests DISCOVERY_MODE PRE_TEST)
endif()

if(MDBRIDGE_BUILD_APP)
    # the main window driven by action scripts (tests/data/gui, format in tests/data/README.md "gui"):
    # saved files against MDATK's, screenshots to $MDBRIDGE_GUI_SCREENSHOTS (default
    # <build>/gui-screenshots; CI artifacts, never compared). Qt Test on the offscreen platform.
    add_executable(mdbridge_gui_tests gui/test_gui.cpp gui/gui_harness.cpp gui/gui_harness.hpp gui/page_host.cpp
                   gui/page_host.hpp)
    set_target_properties(mdbridge_gui_tests PROPERTIES AUTOMOC ON)
    target_link_libraries(mdbridge_gui_tests PRIVATE mdbridge_app mdbridge::render_panels_vulkan Qt6::Test)
    target_compile_definitions(mdbridge_gui_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                               MDBRIDGE_GUI_SCREENSHOT_DIR="${CMAKE_BINARY_DIR}/gui-screenshots")
    mdbridge_target_warnings(mdbridge_gui_tests)
    add_test(NAME gui.cases COMMAND mdbridge_gui_tests)
    set_tests_properties(gui.cases PROPERTIES LABELS gui ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    # the settings pages (src/app/gui/pages, milestone 5 W1-B): MDATK's sections and fields, the group
    # drop-downs, group renames, Groups from chains, reveal; a screenshot of each page (artifacts)
    add_executable(mdbridge_gui_pages_tests gui/test_pages.cpp gui/gui_harness.cpp gui/gui_harness.hpp gui/page_host.cpp
                   gui/page_host.hpp)
    set_target_properties(mdbridge_gui_pages_tests PROPERTIES AUTOMOC ON)
    target_link_libraries(mdbridge_gui_pages_tests PRIVATE mdbridge_app mdbridge::render_panels_vulkan Qt6::Test)
    target_compile_definitions(mdbridge_gui_pages_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
                               MDBRIDGE_GUI_SCREENSHOT_DIR="${CMAKE_BINARY_DIR}/gui-screenshots")
    mdbridge_target_warnings(mdbridge_gui_pages_tests)
    add_test(NAME gui.pages COMMAND mdbridge_gui_pages_tests)
    set_tests_properties(gui.pages PROPERTIES LABELS gui ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    # the GUI's bound field widgets (src/app/gui/forms, milestone 5 W0-C): every field of the fixture
    # projects loaded and saved without an edit byte for byte, edits, undo / redo, AP-01 .. AP-05,
    # AP-15, the mouse wheel. Qt Test on the offscreen platform.
    add_executable(mdbridge_gui_forms_tests gui/test_forms.cpp)
    set_target_properties(mdbridge_gui_forms_tests PROPERTIES AUTOMOC ON)
    target_include_directories(mdbridge_gui_forms_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
    target_link_libraries(mdbridge_gui_forms_tests PRIVATE mdbridge_app Qt6::Test)
    target_compile_definitions(mdbridge_gui_forms_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
    mdbridge_target_warnings(mdbridge_gui_forms_tests)
    add_test(NAME gui.forms COMMAND mdbridge_gui_forms_tests)
    set_tests_properties(gui.forms PROPERTIES LABELS gui ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    # the problems presentation (src/app/gui/problems, milestone 5 W1-D): every validate case's problem
    # on its field of a factory-built form with its location, page badges, jumps, checking as you type
    add_executable(mdbridge_gui_problems_tests gui/test_problems.cpp)
    set_target_properties(mdbridge_gui_problems_tests PROPERTIES AUTOMOC ON)
    target_include_directories(mdbridge_gui_problems_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/core")
    target_link_libraries(mdbridge_gui_problems_tests PRIVATE mdbridge_app Qt6::Test)
    target_compile_definitions(mdbridge_gui_problems_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
    mdbridge_target_warnings(mdbridge_gui_problems_tests)
    add_test(NAME gui.problems COMMAND mdbridge_gui_problems_tests)
    set_tests_properties(gui.problems PROPERTIES LABELS gui ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

    # the interface's job runner and log view (src/app/gui/jobs; Qt Test, offscreen platform): job
    # texts, one job at a time, progress throttle, cancel paths, failures, close-while-running
    # (AP-14), a bounded batched log that keeps the event loop responsive under 10^6 lines
    add_executable(mdbridge_gui_jobs_tests gui/test_jobs.cpp)
    set_target_properties(mdbridge_gui_jobs_tests PROPERTIES AUTOMOC ON)
    target_link_libraries(mdbridge_gui_jobs_tests PRIVATE mdbridge_app Qt6::Test)
    mdbridge_target_warnings(mdbridge_gui_jobs_tests)
    add_test(NAME gui.jobs COMMAND mdbridge_gui_jobs_tests)
    set_tests_properties(gui.jobs PROPERTIES LABELS gui ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 600)

    # every GUI test program links mdbridge_app, which links the video encoder: the FFmpeg DLLs next to it
    if(WIN32 AND MDBRIDGE_VIDEO_RUNTIME_FILES)
        foreach(gui_test mdbridge_gui_tests mdbridge_gui_pages_tests mdbridge_gui_forms_tests mdbridge_gui_problems_tests mdbridge_gui_jobs_tests)
            add_custom_command(TARGET ${gui_test} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MDBRIDGE_VIDEO_RUNTIME_FILES} $<TARGET_FILE_DIR:${gui_test}>
                VERBATIM)
        endforeach()
    endif()
endif()
