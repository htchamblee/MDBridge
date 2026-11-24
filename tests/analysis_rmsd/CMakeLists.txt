# analysis: the optimal-superposition RMSD against MDAnalysis (data/analysis/rmsd_cases.json); its
# own executable until it is folded into the analysis tests
add_executable(mdbridge_analysis_rmsd_tests test_rmsd_qcp.cpp)
target_include_directories(mdbridge_analysis_rmsd_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../core")
target_link_libraries(mdbridge_analysis_rmsd_tests PRIVATE mdbridge::analysis_rmsd mdbridge::io GTest::gtest_main)
target_compile_definitions(mdbridge_analysis_rmsd_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data")
mdbridge_target_warnings(mdbridge_analysis_rmsd_tests)
gtest_discover_tests(mdbridge_analysis_rmsd_tests DISCOVERY_MODE PRE_TEST)
