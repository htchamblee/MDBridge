# render: the 3D structure view (Vulkan 1.3 with a software fallback). Qt GUI, no widgets.
if(NOT QT_FEATURE_vulkan)
    message(FATAL_ERROR "Qt ${Qt6_VERSION} at ${Qt6_DIR} was built without Vulkan support; MDBridge needs a Qt "
                        "with Vulkan (the official Qt builds and most distributions' Qt have it)")
endif()
if(NOT TARGET Vulkan::Headers)  # Qt's package defines it when it found installed Vulkan headers
    FetchContent_MakeAvailable(VulkanHeaders)
endif()

add_library(mdbridge_render STATIC
    include/mdbridge/render/view.hpp
    include/mdbridge/render/vulkan_probe.hpp
    scene.hpp
    shaders/spirv.hpp
    software_view.cpp
    software_view.hpp
    view.cpp
    vulkan_probe.cpp
    vulkan_view.cpp
    vulkan_view.hpp)
add_library(mdbridge::render ALIAS mdbridge_render)
set_target_properties(mdbridge_render PROPERTIES AUTOMOC ON)
target_include_directories(mdbridge_render
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
# Vulkan headers only: Qt and the probe load the system's Vulkan loader at run time
target_link_libraries(mdbridge_render PUBLIC Qt6::Gui PRIVATE Vulkan::Headers)
mdbridge_target_warnings(mdbridge_render)
