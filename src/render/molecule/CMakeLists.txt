# render/molecule: the 3D representations' geometry (MDATK's mdatk/render/molecule/topology.py,
# geometry.py, surface.py, scene.py), the built-in software renderer (software.py), the density shadow (mdatk/density.py: project_density,
# shadow_rgba, anchored) and the pieces of scipy.ndimage they use. No Qt: used by the
# software renderer, the Vulkan renderer and the structure panel, and tested without a GUI.
add_library(mdbridge_render_molecule STATIC
    density_shadow.cpp
    geometry.cpp
    include/mdbridge/render/molecule/density_shadow.hpp
    include/mdbridge/render/molecule/geometry.hpp
    include/mdbridge/render/molecule/ndimage.hpp
    include/mdbridge/render/molecule/scene.hpp
    include/mdbridge/render/molecule/software.hpp
    include/mdbridge/render/molecule/surface.hpp
    include/mdbridge/render/molecule/topology.hpp
    ndimage.cpp
    scene.cpp
    software.cpp
    surface.cpp
    topology.cpp)
add_library(mdbridge::render_molecule ALIAS mdbridge_render_molecule)
target_include_directories(mdbridge_render_molecule PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
# io: the density map (DensityMap) in the shadow's interface; analysis: kabsch for anchored()
find_package(Threads REQUIRED)  # the software renderer's worker threads (std::thread)
target_link_libraries(mdbridge_render_molecule PUBLIC mdbridge::core mdbridge::io PRIVATE mdbridge::analysis
                      Threads::Threads)
target_compile_features(mdbridge_render_molecule PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_render_molecule)
# the surface's float32 grid and the software renderer's float32 pixels are rounded as numpy / scipy
# round them: no fused multiply-adds
if(MSVC)
    target_compile_options(mdbridge_render_molecule PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_render_molecule PRIVATE -ffp-contract=off)
endif()
