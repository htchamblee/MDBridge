# chem: nucleic-acid chemistry ported from MDATK: PDB hygiene (pdbtools.cpp), secondary structure
# (dot-bracket notation and the builtin base-pair detection), base-pair annotation (annotate.cpp),
# SMILES, molecular graphs and the restraint dictionaries (smiles.cpp, dictionaries.cpp), the
# structure builder's templates, sequence functions and Watson-Crick helices (helix.cpp), the
# rigid-body density fit (density_fit.cpp, on src/numerics) and the generated reference tables of
# the builder and the restraint dictionaries (reference_tables.cpp: helix_tables.inc,
# component_tables.inc from tools/gen_chem_tables.py). Residue names come from io's
# parseResidueName, mmCIF from io's cif.
add_library(mdbridge_chem STATIC
    annotate.cpp
    density_fit.cpp
    dictionaries.cpp
    helix.cpp
    pdb_names.cpp
    pdbtools.cpp
    reference_tables.cpp
    secstruct.cpp
    smiles.cpp)
add_library(mdbridge::chem ALIAS mdbridge_chem)
target_include_directories(mdbridge_chem PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_chem PUBLIC mdbridge::core mdbridge::io PRIVATE mdbridge::numerics)
target_compile_features(mdbridge_chem PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_chem)
# distances (secstruct, pdbtools), chiral volumes / torsions (dictionaries) and helix coordinates
# (helix) reproduced bit for bit (numpy's and Python's order of operations), and the density fit's
# numbers the same with every compiler: no contraction into fused multiply-adds
if(MSVC)
    set_source_files_properties(secstruct.cpp pdbtools.cpp dictionaries.cpp helix.cpp density_fit.cpp
        PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(secstruct.cpp pdbtools.cpp dictionaries.cpp helix.cpp density_fit.cpp
        PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
