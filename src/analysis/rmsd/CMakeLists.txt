# analysis (RMSD only, for now): the optimal-superposition RMSD MDATK's rmsd analysis gets from
# MDAnalysis (docs/specs/m3-rmsd.md). A target of its own until it is folded into the analysis
# library.
set(MDBRIDGE_RMSD_QCP_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../rmsd_qcp.cpp")
add_library(mdbridge_analysis_rmsd STATIC "${MDBRIDGE_RMSD_QCP_SOURCE}")
add_library(mdbridge::analysis_rmsd ALIAS mdbridge_analysis_rmsd)
target_include_directories(mdbridge_analysis_rmsd PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(mdbridge_analysis_rmsd PUBLIC mdbridge::core)
target_compile_features(mdbridge_analysis_rmsd PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_analysis_rmsd)
# MDAnalysis' results reproduced bit for bit (also the chaotic degenerate cases): every operation
# rounded, no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties("${MDBRIDGE_RMSD_QCP_SOURCE}" PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties("${MDBRIDGE_RMSD_QCP_SOURCE}" PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
