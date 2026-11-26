# render/panels: the panels of MDATK's multi-panel figure (mdatk/render/panels.py, seqpanels.py) on
# the figure model of render core, and build_figure (mdatk/render/video.py). No Qt: shared by the
# video loop, the charts, the preview and the GUI. One source per panel family (panel.hpp); the
# standalone charts (mdatk/render/charts.py) in charts.cpp, written by a FigureWriter (Qt: figure_qt.hpp).
add_library(mdbridge_render_panels STATIC
    include/mdbridge/render/panels/build.hpp
    include/mdbridge/render/panels/panel.hpp
    build.cpp
    panel.cpp
    registry.cpp
    include/mdbridge/render/panels/plots.hpp
    plots.cpp
    include/mdbridge/render/panels/seqpanels.hpp
    seqpanels.cpp
    include/mdbridge/render/panels/structure.hpp
    structure.cpp
    include/mdbridge/render/panels/charts.hpp
    charts.cpp
    include/mdbridge/render/panels/checks.hpp
    checks.cpp)
add_library(mdbridge::render_panels ALIAS mdbridge_render_panels)
target_include_directories(mdbridge_render_panels PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_render_panels PUBLIC mdbridge::render_core mdbridge::render_molecule mdbridge::io mdbridge::core
    PRIVATE mdbridge::chem mdbridge::numerics)  # chem: dot-bracket (secondary panel); io: density maps (structure
                                                # view); numerics: density_centroid
target_compile_features(mdbridge_render_panels PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_panels)
# the residue panels' and the structure view's arithmetic rounds as numpy rounds it: no contraction
# into fused multiply-adds
if(MSVC)
    set_source_files_properties(seqpanels.cpp structure.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(seqpanels.cpp structure.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()

# the structure view's 3D layer on the GPU (Vulkan scene renderer) and the renderer choice that
# installs it (vulkan_stage.hpp); a library of its own so the panels need no Vulkan
add_library(mdbridge_render_panels_vulkan STATIC
    include/mdbridge/render/panels/vulkan_stage.hpp
    vulkan_stage.cpp)
add_library(mdbridge::render_panels_vulkan ALIAS mdbridge_render_panels_vulkan)
target_link_libraries(mdbridge_render_panels_vulkan PUBLIC mdbridge::render_panels mdbridge::render_vulkan_scene)
mdbridge_target_warnings(mdbridge_render_panels_vulkan)
