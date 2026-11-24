# analysis spatial helpers (W0.7): MDATK's mdatk/spatial.py. A target of its own until the
# analysis library (src/analysis/CMakeLists.txt) takes spatial.cpp in.
set(MDBRIDGE_ANALYSIS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_library(mdbridge_analysis_spatial STATIC "${MDBRIDGE_ANALYSIS_DIR}/spatial.cpp")
add_library(mdbridge::analysis_spatial ALIAS mdbridge_analysis_spatial)
target_include_directories(mdbridge_analysis_spatial PUBLIC "${MDBRIDGE_ANALYSIS_DIR}/include")
target_link_libraries(mdbridge_analysis_spatial PUBLIC mdbridge::core)
target_compile_features(mdbridge_analysis_spatial PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_analysis_spatial)
# distances and products rounded as scipy / numpy round them: no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties("${MDBRIDGE_ANALYSIS_DIR}/spatial.cpp" PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties("${MDBRIDGE_ANALYSIS_DIR}/spatial.cpp" PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
