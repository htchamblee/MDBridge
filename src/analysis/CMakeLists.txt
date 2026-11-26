# analysis: MDATK's analyses (mdatk/analysis/*.py): the framework (registry, instance ids, cache
# keys, context, single-pass frame-parallel runner), the analyses themselves, and the analysis cache and analyze
# pipeline (mdatk/io/cache.py, mdatk/pipeline.py; .npz entries written with zlib).
add_library(mdbridge_analysis STATIC
    analysis.cpp
    basepairs.cpp
    builtin.hpp
    cache.cpp
    cleaning.cpp
    common.cpp
    contacts.cpp
    context.cpp
    energy.cpp
    energy_parameters.cpp
    hbonds.cpp
    keys.cpp
    mg.cpp
    pipeline.cpp
    registry.cpp
    rmsd.cpp
    rmsd_qcp.cpp
    runner.cpp
    spatial.cpp
    stacking.cpp
    structure.cpp)
add_library(mdbridge::analysis ALIAS mdbridge_analysis)
find_package(Threads REQUIRED)  # the runner's worker threads (std::thread)
target_include_directories(mdbridge_analysis PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_analysis PUBLIC mdbridge::core mdbridge::io PRIVATE mdbridge::chem Threads::Threads ${MDBRIDGE_ZLIB_TARGET})
target_compile_features(mdbridge_analysis PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_analysis)
# numerics reproduced as numpy rounds them: no contraction into fused multiply-adds
if(MSVC)
    target_compile_options(mdbridge_analysis PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_analysis PRIVATE -ffp-contract=off)
endif()
