# app: the Qt interface and the command line, two thin front ends over the libraries.
FetchContent_MakeAvailable(CLI11)

add_library(mdbridge_app STATIC
    cli.cpp
    cli.hpp
    gui/gui.cpp
    gui/gui.hpp
    gui/main_window.cpp
    gui/main_window.hpp
    gui/project_panel.cpp
    gui/project_panel.hpp
    project_summary.cpp
    project_summary.hpp)
set_target_properties(mdbridge_app PROPERTIES AUTOMOC ON)
target_include_directories(mdbridge_app PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_app
    PUBLIC mdbridge::core mdbridge::render Qt6::Widgets
    PRIVATE CLI11::CLI11)
mdbridge_target_warnings(mdbridge_app)

set(MDBRIDGE_ICON_ICO "${PROJECT_SOURCE_DIR}/data/icons/mdbridge.ico")
configure_file(windows/mdbridge.rc.in "${CMAKE_CURRENT_BINARY_DIR}/mdbridge.rc" @ONLY)

# mdbridge_program(<name> [GUI]): a front end; GUI = Windows program without a console window
function(mdbridge_program name)
    cmake_parse_arguments(ARG "GUI" "" "" ${ARGN})
    add_executable(${name} main.cpp)
    target_link_libraries(${name} PRIVATE mdbridge_app)
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
