# numerics: the pieces of SciPy that MDATK's density fit, density shadow, molecular surface and helix
# builder use, ported or vendored with their notices (LICENSING.md): scipy.ndimage filters and
# interpolation (ndimage.cpp), scipy.fft.rfftn / irfftn on SciPy's own FFT (fft.cpp, third_party/ducc0),
# scipy.optimize's Powell minimiser with Brent's line search (powell.cpp) and least_squares' trust
# region reflective method (least_squares.cpp, lsmr.cpp), cKDTree.query_pairs (pairs.cpp), and
# mdatk/density.py's density_centroid (density.cpp). No Qt; chem, analysis and render may use it.

# ducc0 0.40.0's FFT as SciPy 1.18.1 vendors it (tools/vendor_ducc0.py): third-party code, compiled
# with the compiler's default warnings and never edited; every file is checked against SHA256SUMS
set(MDBRIDGE_DUCC0_DIR "${PROJECT_SOURCE_DIR}/third_party/ducc0")
file(STRINGS "${MDBRIDGE_DUCC0_DIR}/SHA256SUMS" ducc0_sums)
foreach(line IN LISTS ducc0_sums)
    string(REGEX MATCH "^([0-9a-f]+)  (.+)$" matched "${line}")
    file(SHA256 "${MDBRIDGE_DUCC0_DIR}/${CMAKE_MATCH_2}" actual)
    if(NOT actual STREQUAL CMAKE_MATCH_1)
        message(FATAL_ERROR "third_party/ducc0/${CMAKE_MATCH_2} differs from SciPy 1.18.1's copy "
                            "(SHA256SUMS; recreate it with tools/vendor_ducc0.py)")
    endif()
endforeach()
add_library(mdbridge_ducc0 STATIC
    "${MDBRIDGE_DUCC0_DIR}/ducc0/infra/mav.cc"
    "${MDBRIDGE_DUCC0_DIR}/ducc0/infra/string_utils.cc"
    "${MDBRIDGE_DUCC0_DIR}/ducc0/infra/threading.cc")
target_include_directories(mdbridge_ducc0 SYSTEM PUBLIC "${MDBRIDGE_DUCC0_DIR}")
find_package(Threads REQUIRED)  # ducc0's thread pool
target_link_libraries(mdbridge_ducc0 PUBLIC Threads::Threads)
target_compile_features(mdbridge_ducc0 PUBLIC cxx_std_17)
set_target_properties(mdbridge_ducc0 PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MINGW)
    # as SciPy builds it: no low-level threading with mingw-w64 (no pthread_atfork)
    target_compile_definitions(mdbridge_ducc0 PUBLIC DUCC0_NO_LOWLEVEL_THREADING)
endif()
if(MSVC)
    target_compile_options(mdbridge_ducc0 PRIVATE /W0 /fp:precise)
else()
    # SciPy: -Wno-psabi (notes about the GCC 10.1 ABI of SIMD types passed by value)
    target_compile_options(mdbridge_ducc0 PUBLIC $<$<CXX_COMPILER_ID:GNU>:-Wno-psabi> PRIVATE -w -ffp-contract=off)
endif()

add_library(mdbridge_numerics STATIC
    density.cpp
    fft.cpp
    include/mdbridge/numerics/density.hpp
    include/mdbridge/numerics/fft.hpp
    include/mdbridge/numerics/least_squares.hpp
    include/mdbridge/numerics/ndimage.hpp
    include/mdbridge/numerics/pairs.hpp
    include/mdbridge/numerics/powell.hpp
    least_squares.cpp
    lsmr.cpp
    ndimage.cpp
    numerics_internal.hpp
    pairs.cpp
    powell.cpp)
add_library(mdbridge::numerics ALIAS mdbridge_numerics)
target_include_directories(mdbridge_numerics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
# io: DensityMap (density_centroid); core: numpy's pairwise sum, the cell grid of the pair search
find_package(Threads REQUIRED)  # the finite differences' optional worker threads
target_link_libraries(mdbridge_numerics PUBLIC mdbridge::core mdbridge::io PRIVATE mdbridge_ducc0 Threads::Threads)
target_compile_features(mdbridge_numerics PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_numerics)
# every product and sum rounds as SciPy's: no contraction into fused multiply-adds
if(MSVC)
    target_compile_options(mdbridge_numerics PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_numerics PRIVATE -ffp-contract=off)
endif()
