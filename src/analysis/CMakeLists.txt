# analysis: MDATK's analyses (mdatk/analysis/*.py): the framework (registry, instance ids, cache
# keys, context, single-pass runner) and the analyses themselves.
add_library(mdbridge_analysis STATIC
    analysis.cpp
    basepairs.cpp
    builtin.hpp
    common.cpp
    contacts.cpp
    context.cpp
    energy.cpp
    hbonds.cpp
    keys.cpp
    mg.cpp
    registry.cpp
    rmsd.cpp
    runner.cpp
    spatial.cpp
    stacking.cpp
    structure.cpp)
add_library(mdbridge::analysis ALIAS mdbridge_analysis)
target_include_directories(mdbridge_analysis PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_analysis PUBLIC mdbridge::core mdbridge::io)
target_compile_features(mdbridge_analysis PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_analysis)
# numerics reproduced as numpy rounds them: no contraction into fused multiply-adds
if(MSVC)
    target_compile_options(mdbridge_analysis PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_analysis PRIVATE -ffp-contract=off)
endif()
