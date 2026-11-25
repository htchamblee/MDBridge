# render core: the Qt-free part of the render module (colours, colormaps, themes, ticks, fonts and
# text metrics, residue shading, orientation, the example dataset with numpy's PCG64), testable
# without Qt. Namespace mdbridge::render.
add_library(mdbridge_render_core STATIC
    include/mdbridge/render/color.hpp
    include/mdbridge/render/colormap.hpp
    include/mdbridge/render/font.hpp
    include/mdbridge/render/text.hpp
    include/mdbridge/render/theme.hpp
    color.cpp
    colormap.cpp
    font.cpp
    mpl_colors.inc
    mpl_tables.hpp
    include/mdbridge/render/ticks.hpp
    text.cpp
    theme.cpp
    ticks.cpp
    include/mdbridge/render/example.hpp
    include/mdbridge/render/numpy_random.hpp
    include/mdbridge/render/orientation.hpp
    include/mdbridge/render/shading.hpp
    example.cpp
    numpy_random.cpp
    orientation.cpp
    shading.cpp
    ziggurat_constants.inc)
add_library(mdbridge::render_core ALIAS mdbridge_render_core)
target_include_directories(mdbridge_render_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mdbridge_render_core PUBLIC mdbridge::core PRIVATE mdbridge::analysis mdbridge::chem mdbridge::io)
target_compile_features(mdbridge_render_core PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_core)
# lookup tables, linspace and tick arithmetic rounded as numpy / Python round them: no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties(colormap.cpp ticks.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
    set_source_files_properties(example.cpp numpy_random.cpp orientation.cpp shading.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(colormap.cpp ticks.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
    set_source_files_properties(example.cpp numpy_random.cpp orientation.cpp shading.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
