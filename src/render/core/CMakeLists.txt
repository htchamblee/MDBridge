# render core: the Qt-free part of the render module (colours, colormaps, themes, ticks; later
# molecule geometry, shading ...), testable without Qt. Namespace mdbridge::render.
add_library(mdbridge_render_core STATIC
    include/mdbridge/render/color.hpp
    include/mdbridge/render/colormap.hpp
    include/mdbridge/render/theme.hpp
    color.cpp
    colormap.cpp
    mpl_colors.inc
    mpl_tables.hpp
    include/mdbridge/render/ticks.hpp
    theme.cpp
    ticks.cpp)
add_library(mdbridge::render_core ALIAS mdbridge_render_core)
target_include_directories(mdbridge_render_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_render_core PUBLIC mdbridge::core)
target_compile_features(mdbridge_render_core PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_core)
# lookup tables, linspace and tick arithmetic rounded as numpy / Python round them: no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties(colormap.cpp ticks.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(colormap.cpp ticks.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
