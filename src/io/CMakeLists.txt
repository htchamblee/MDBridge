# io: reading and writing the files MDATK reads and writes (topologies, trajectories, PDB, maps,
# restraint dictionaries, external tables), with MDAnalysis' conventions where MDATK relies on them.
add_library(mdbridge_io STATIC
    amber_ascii.cpp
    dcd.cpp
    external.cpp
    gromacs.cpp
    loader.cpp
    mrc.cpp
    netcdf.cpp
    parm7.cpp
    pdb.cpp
    pdb_export.cpp
    pdb_topology.cpp
    readers.hpp
    trajectory.cpp)
add_library(mdbridge::io ALIAS mdbridge_io)
target_include_directories(mdbridge_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mdbridge_io PUBLIC mdbridge::core PRIVATE ${MDBRIDGE_ZLIB_TARGET})
target_compile_features(mdbridge_io PUBLIC cxx_std_20)
mdbridge_target_warnings(mdbridge_io)
# numerics reproduced bit for bit (numpy's order of operations): no contraction into fused
# multiply-adds
if(MSVC)
    set_source_files_properties(mrc.cpp PROPERTIES COMPILE_OPTIONS /fp:precise)
else()
    set_source_files_properties(mrc.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
endif()
