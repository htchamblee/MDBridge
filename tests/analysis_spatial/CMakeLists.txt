# analysis spatial helpers against MDATK (data/analysis/spatial_extra_cases.json) and the periodic
# pair search against MDAnalysis (data/analysis/distance_cases.json)
add_executable(mdbridge_analysis_spatial_tests test_spatial.cpp test_periodic.cpp)
target_include_directories(mdbridge_analysis_spatial_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../core")
target_link_libraries(mdbridge_analysis_spatial_tests PRIVATE mdbridge::analysis GTest::gtest_main)
target_compile_definitions(mdbridge_analysis_spatial_tests PRIVATE MDBRIDGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data")
mdbridge_target_warnings(mdbridge_analysis_spatial_tests)
# the expected sums and products computed here as numpy does: no fused multiply-adds
if(MSVC)
    target_compile_options(mdbridge_analysis_spatial_tests PRIVATE /fp:precise)
else()
    target_compile_options(mdbridge_analysis_spatial_tests PRIVATE -ffp-contract=off)
endif()
gtest_discover_tests(mdbridge_analysis_spatial_tests DISCOVERY_MODE PRE_TEST)
