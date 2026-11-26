# chem: nucleic-acid chemistry ported from MDATK (secondary structure: dot-bracket notation and the
# builtin base-pair detection; SMILES, molecular graphs and the restraint dictionaries) and the
# generated reference tables of the builder and the restraint dictionaries (reference_tables.cpp:
# helix_tables.inc, component_tables.inc from tools/gen_chem_tables.py). Residue names come from io's
# parseResidueName, mmCIF from io's cif.
add_library(mdbridge_chem STATIC
    dictionaries.cpp
    reference_tables.cpp
    secstruct.cpp
    smiles.cpp)
add_library(mdbridge::chem ALIAS mdbridge_chem)
target_include_directories(mdbridge_chem PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_chem PUBLIC mdbridge::core mdbridge::io)
target_compile_features(mdbridge_chem PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_chem)
# distances (secstruct) and chiral volumes / torsions (dictionaries) reproduced bit for bit (numpy's
# and Python's order of operations): no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties(secstruct.cpp dictionaries.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(secstruct.cpp dictionaries.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
