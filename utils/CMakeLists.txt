project(utils)

add_definitions(-DTEST)
add_executable(sp_topo30_mesh sp_srtm_rel.cpp)

if (OPENMP_FOUND)
set_target_properties(sp_topo30_mesh PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif ()
