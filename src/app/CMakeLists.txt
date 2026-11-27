# app: the Qt interface and the command line, two thin front ends over the libraries. The command
# line's MDATK commands are in commands/ (mdbridge::commands, no Qt).

add_library(mdbridge_app STATIC
    chain_groups.cpp
    chain_groups.hpp
    cli.cpp
    cli.hpp
    fonts_command.cpp
    fonts_command.hpp
    loaded_libraries.cpp
    loaded_libraries.hpp
    gui/actions/window_actions.cpp
    gui/actions/window_actions.hpp
    gui/background.hpp
    gui/dialogs/image_dialog.cpp
    gui/dialogs/image_dialog.hpp
    gui/dialogs/prompts.cpp
    gui/dialogs/prompts.hpp
    gui/dialogs/text_dialog.cpp
    gui/dialogs/text_dialog.hpp
    gui/display_path.cpp
    gui/display_path.hpp
    gui/forms/custom_fields.cpp
    gui/forms/custom_fields.hpp
    gui/forms/field_factory.cpp
    gui/forms/field_factory.hpp
    gui/forms/field_spec.cpp
    gui/forms/field_spec.hpp
    gui/forms/field_widgets.cpp
    gui/forms/field_widgets.hpp
    gui/forms/flow_layout.cpp
    gui/forms/flow_layout.hpp
    gui/forms/form_context.cpp
    gui/forms/form_context.hpp
    gui/forms/list_editor.cpp
    gui/forms/list_editor.hpp
    gui/gui.cpp
    gui/gui.hpp
    gui/jobs/job_runner.cpp
    gui/jobs/job_runner.hpp
    gui/jobs/log_view.cpp
    gui/jobs/log_view.hpp
    gui/log_capture.hpp
    gui/main_window.cpp
    gui/main_window.hpp
    gui/pages/analyses_page.cpp
    gui/pages/analyses_page.hpp
    gui/pages/form_page.cpp
    gui/pages/form_page.hpp
    gui/pages/groups_page.cpp
    gui/pages/groups_page.hpp
    gui/pages/input_page.cpp
    gui/pages/input_page.hpp
    gui/pages/layout_page.cpp
    gui/pages/layout_page.hpp
    gui/pages/outputs_page.cpp
    gui/pages/outputs_page.hpp
    gui/pages/page_registry.cpp
    gui/pages/page_registry.hpp
    gui/pages/placeholder_page.cpp
    gui/pages/placeholder_page.hpp
    gui/pages/structure_page.cpp
    gui/pages/structure_page.hpp
    gui/preview_widget.cpp
    gui/preview_widget.hpp
    gui/problems/problems_panel.cpp
    gui/problems/problems_panel.hpp
    gui/script_runner.cpp
    gui/script_runner.hpp
    gui/tools/builder_window.cpp
    gui/tools/builder_window.hpp
    gui/tools/host.hpp
    gui/tools/pdb_tools_window.cpp
    gui/tools/pdb_tools_window.hpp
    gui/tools/tool_actions.cpp
    gui/tools/tool_actions.hpp
    gui/tools/tool_window.cpp
    gui/tools/tool_window.hpp
    gui/tools/tool_windows.cpp
    gui/tools/tool_windows.hpp
    project_summary.cpp
    project_summary.hpp
    pubchem_qt.cpp
    pubchem_qt.hpp
    render_backend_qt.cpp
    render_backend_qt.hpp)
set_target_properties(mdbridge_app PROPERTIES AUTOMOC ON)
target_include_directories(mdbridge_app PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_app
    PUBLIC mdbridge::core mdbridge::render Qt6::Widgets
    PRIVATE Qt6::Network mdbridge::analysis mdbridge::chem mdbridge::commands mdbridge::io mdbridge::numerics mdbridge::render_panels mdbridge::render_panels_vulkan
            mdbridge::video CLI11::CLI11)
mdbridge_target_warnings(mdbridge_app)

set(MDBRIDGE_ICON_ICO "${PROJECT_SOURCE_DIR}/data/icons/mdbridge.ico")
configure_file(windows/mdbridge.rc.in "${CMAKE_CURRENT_BINARY_DIR}/mdbridge.rc" @ONLY)

# mdbridge_program(<name> [GUI]): a front end; GUI = Windows program without a console window
function(mdbridge_program name)
    cmake_parse_arguments(ARG "GUI" "" "" ${ARGN})
    add_executable(${name} main.cpp)
    target_link_libraries(${name} PRIVATE mdbridge_app mdbridge::render_panels_vulkan)
    qt_add_resources(${name} icons PREFIX /mdbridge BASE "${PROJECT_SOURCE_DIR}/data"
                     FILES "${PROJECT_SOURCE_DIR}/data/icons/mdbridge.png")
    if(WIN32)
        target_sources(${name} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/mdbridge.rc")
    endif()
    if(ARG_GUI)
        target_compile_definitions(${name} PRIVATE MDBRIDGE_GUI_ENTRY)
        set_target_properties(${name} PROPERTIES WIN32_EXECUTABLE ON)
    endif()
    mdbridge_target_warnings(${name})
    if(WIN32 AND MDBRIDGE_VIDEO_RUNTIME_FILES)  # the FFmpeg DLLs next to the build tree's program
        add_custom_command(TARGET ${name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MDBRIDGE_VIDEO_RUNTIME_FILES} $<TARGET_FILE_DIR:${name}>
            VERBATIM)
    endif()
    install(TARGETS ${name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endfunction()

mdbridge_program(mdbridge)
if(WIN32)
    mdbridge_program(mdbridge-gui GUI)
endif()

if(WIN32)
    # Qt's DLLs and plugins next to the programs (windeployqt), and the MSVC runtime app-local, so
    # the installer needs no administrator rights and no separate runtime installer
    qt_generate_deploy_app_script(TARGET mdbridge-gui OUTPUT_SCRIPT deploy_script NO_UNSUPPORTED_PLATFORM_ERROR
                                  DEPLOY_TOOL_OPTIONS --no-opengl-sw --no-system-d3d-compiler --no-translations)
    install(SCRIPT "${deploy_script}")
    set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
    include(InstallRequiredSystemLibraries)
else()
    install(FILES "${PROJECT_SOURCE_DIR}/data/mdbridge.desktop" DESTINATION ${CMAKE_INSTALL_DATADIR}/applications)
    install(FILES "${PROJECT_SOURCE_DIR}/data/icons/mdbridge.png"
            DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/256x256/apps)
endif()
