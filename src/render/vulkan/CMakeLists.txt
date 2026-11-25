# render/vulkan: Vulkan 1.3 without a window (device and queue, offscreen render target and
# read-back, pipelines and the pipeline cache) and the Vulkan probe. No Qt: the batch renderer and
# the tests use it without a GUI; the GUI's window keeps QVulkanWindow (src/render). The system's
# Vulkan loader is opened at run time (never linked or bundled); only the headers are needed.
if(NOT TARGET Vulkan::Headers)
    FetchContent_MakeAvailable(VulkanHeaders)
endif()

add_library(mdbridge_render_vulkan STATIC
    context.cpp
    include/mdbridge/render/vulkan/context.hpp
    include/mdbridge/render/vulkan/functions.hpp
    include/mdbridge/render/vulkan/loader.hpp
    include/mdbridge/render/vulkan/offscreen.hpp
    include/mdbridge/render/vulkan/pipeline.hpp
    include/mdbridge/render/vulkan_probe.hpp
    loader.cpp
    offscreen.cpp
    pipeline.cpp
    vulkan_probe.cpp)
add_library(mdbridge::render_vulkan ALIAS mdbridge_render_vulkan)
target_include_directories(mdbridge_render_vulkan PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_render_vulkan PUBLIC Vulkan::Headers PRIVATE mdbridge::core ${CMAKE_DL_LIBS})
target_compile_definitions(mdbridge_render_vulkan PUBLIC VK_NO_PROTOTYPES)
target_compile_features(mdbridge_render_vulkan PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_vulkan)

# The GPU renderer of the molecule scene (impostor spheres and sticks, meshes, fog, glow, outlines,
# translucency, anti-aliasing) on the base above: offscreen for batch work, into a window's command
# buffer for the GUI. No Qt. Shaders: src/render/shaders/molecule_* (tools/compile_shaders.py).
add_library(mdbridge_render_vulkan_scene STATIC
    include/mdbridge/render/vulkan/scene_renderer.hpp
    scene_renderer.cpp)
add_library(mdbridge::render_vulkan_scene ALIAS mdbridge_render_vulkan_scene)
target_include_directories(mdbridge_render_vulkan_scene PRIVATE "${PROJECT_SOURCE_DIR}/src/render")  # shaders/spirv.hpp
target_link_libraries(mdbridge_render_vulkan_scene PUBLIC mdbridge::render_vulkan mdbridge::render_molecule)
mdbridge_target_warnings(mdbridge_render_vulkan_scene)
# primitives are set up in the software renderer's float32 / double arithmetic: no contraction
if(MSVC)
    target_compile_options(mdbridge_render_vulkan_scene PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_render_vulkan_scene PRIVATE -ffp-contract=off)
endif()
