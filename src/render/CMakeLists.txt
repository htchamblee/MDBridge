# render: the 3D structure view (Vulkan 1.3 with a software fallback) and the figure framework's
# QPainter backends (raster, SVG, PDF). Qt GUI, no widgets.
if(NOT QT_FEATURE_vulkan)
    message(FATAL_ERROR "Qt ${Qt6_VERSION} at ${Qt6_DIR} was built without Vulkan support; MDBridge needs a Qt "
                        "with Vulkan (the official Qt builds and most distributions' Qt have it)")
endif()
add_library(mdbridge_render STATIC
    include/mdbridge/render/figure_qt.hpp
    include/mdbridge/render/text_qt.hpp
    figure_qt.cpp
    include/mdbridge/render/view.hpp
    scene.hpp
    shaders/spirv.hpp
    software_view.cpp
    software_view.hpp
    text_qt.cpp
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
target_link_libraries(mdbridge_render PUBLIC Qt6::Gui Qt6::Svg mdbridge::render_core mdbridge::render_vulkan)
mdbridge_target_warnings(mdbridge_render)

# the figure fonts (data/fonts: DejaVu 2.35 as Matplotlib ships it; licence in packaging/licenses)
# compiled into the programs: figures never depend on the system's fonts
set(MDBRIDGE_FONT_FILES DejaVuSans.ttf DejaVuSans-Bold.ttf DejaVuSans-Oblique.ttf DejaVuSansMono.ttf
                        DejaVuSansMono-Bold.ttf)
list(TRANSFORM MDBRIDGE_FONT_FILES PREPEND "${PROJECT_SOURCE_DIR}/data/fonts/" OUTPUT_VARIABLE mdbridge_font_paths)
qt_add_resources(mdbridge_render mdbridge_fonts PREFIX /mdbridge/fonts BASE "${PROJECT_SOURCE_DIR}/data/fonts"
                 FILES ${mdbridge_font_paths})
