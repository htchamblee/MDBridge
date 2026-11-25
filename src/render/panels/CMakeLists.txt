# render/panels: the panels of MDATK's multi-panel figure (mdatk/render/panels.py, seqpanels.py) on
# the figure model of render core, and build_figure (mdatk/render/video.py). No Qt: shared by the
# video loop, the charts, the preview and the GUI. One source per panel family (panel.hpp).
add_library(mdbridge_render_panels STATIC
    include/mdbridge/render/panels/build.hpp
    include/mdbridge/render/panels/panel.hpp
    include/mdbridge/render/panels/structure.hpp
    build.cpp
    panel.cpp
    registry.cpp
    plots.cpp
    structure.cpp)
add_library(mdbridge::render_panels ALIAS mdbridge_render_panels)
target_include_directories(mdbridge_render_panels PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_render_panels PUBLIC mdbridge::render_core mdbridge::render_molecule mdbridge::io mdbridge::core)
target_compile_features(mdbridge_render_panels PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_panels)
# the structure view's float32 trace arithmetic and view products follow numpy's rounding: no fused
# multiply-adds
if(MSVC)
    set_source_files_properties(structure.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(structure.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
