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
