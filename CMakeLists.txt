project(sph_ex)
cmake_minimum_required(VERSION 2.8)

find_package(OpenMP)
add_definitions(${OpenMP_C_FLAGS})

add_subdirectory(spherepack_c)
add_subdirectory(src)

