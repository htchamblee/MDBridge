# render: the 3D structure view (Vulkan 1.3 with a software fallback) and the figure framework's
# QPainter backends (raster, SVG, PDF). Qt GUI, no widgets.
if(NOT QT_FEATURE_vulkan)
    message(FATAL_ERROR "Qt ${Qt6_VERSION} at ${Qt6_DIR} was built without Vulkan support; MDBridge needs a Qt "
                        "with Vulkan (the official Qt builds and most distributions' Qt have it)")
endif()
add_library(mdbridge_render STATIC
    include/mdbridge/render/figure_qt.hpp
    include/mdbridge/render/text_qt.hpp
    include/mdbridge/render/video_qt.hpp
    figure_qt.cpp
    include/mdbridge/render/view.hpp
    scene.hpp
    shaders/spirv.hpp
    software_view.cpp
    software_view.hpp
    text_qt.cpp
    video_qt.cpp
    view.cpp
    vulkan_view.cpp
    vulkan_view.hpp)
add_library(mdbridge::render ALIAS mdbridge_render)
set_target_properties(mdbridge_render PROPERTIES AUTOMOC ON)
target_include_directories(mdbridge_render
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
# Vulkan headers only: Qt and the probe (render/vulkan) load the system's Vulkan loader at run time
# Qt Svg (LGPL, QSvgGenerator) for SVG figures; QPdfWriter is part of Qt Gui
target_link_libraries(mdbridge_render PUBLIC Qt6::Gui Qt6::Svg mdbridge::render_core mdbridge::render_vulkan
    mdbridge::video_loop  # the video frame loop composed with QPainter (video_qt.hpp)
    mdbridge::render_vulkan_scene)  # the 3D view: the molecule scene on the GPU
mdbridge_target_warnings(mdbridge_render)

# the bundled fonts (data/fonts: IBM Plex Serif, Plex Math and Plex Mono as IBM's releases have
# them, data/sources/fonts/README.md; licence in packaging/licenses/IBM-Plex-OFL-1.1.txt) compiled
# into the programs: figures and the interface never depend on the system's fonts
set(MDBRIDGE_FONT_FILES IBMPlexSerif-Regular.ttf IBMPlexSerif-Bold.ttf IBMPlexSerif-Italic.ttf
                        IBMPlexSerif-BoldItalic.ttf IBMPlexMath-Regular.ttf IBMPlexMono-Regular.ttf
                        IBMPlexMono-Bold.ttf)
list(TRANSFORM MDBRIDGE_FONT_FILES PREPEND "${PROJECT_SOURCE_DIR}/data/fonts/" OUTPUT_VARIABLE mdbridge_font_paths)
qt_add_resources(mdbridge_render mdbridge_fonts PREFIX /mdbridge/fonts BASE "${PROJECT_SOURCE_DIR}/data/fonts"
                 FILES ${mdbridge_font_paths})
