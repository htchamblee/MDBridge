# chem: nucleic-acid chemistry ported from MDATK: PDB hygiene (pdbtools.cpp), secondary structure
# (dot-bracket notation and the builtin base-pair detection) and the generated reference tables of
# the builder and the restraint dictionaries (reference_tables.cpp: helix_tables.inc,
# component_tables.inc from tools/gen_chem_tables.py). Residue names come from io's
# parseResidueName.
add_library(mdbridge_chem STATIC
    pdbtools.cpp
    reference_tables.cpp
    secstruct.cpp)
add_library(mdbridge::chem ALIAS mdbridge_chem)
target_include_directories(mdbridge_chem PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_chem PUBLIC mdbridge::core mdbridge::io)
target_compile_features(mdbridge_chem PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_chem)
# distances reproduced bit for bit (numpy's order of operations): no contraction into fused
# multiply-adds
if(MSVC)
    set_source_files_properties(secstruct.cpp pdbtools.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(secstruct.cpp pdbtools.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
