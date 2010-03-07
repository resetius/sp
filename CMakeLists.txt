project(sph_ex)
cmake_minimum_required(VERSION 2.8)

find_package(OpenMP)
add_definitions(${OpenMP_C_FLAGS})

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/spherepack_c)
add_subdirectory(spherepack_c)
endif ()

add_subdirectory(src)

