project(sph_ex)
cmake_minimum_required(VERSION 2.8)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CheckedAddSubdir)

find_package(OpenMP)
add_definitions(${OpenMP_C_FLAGS})

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/spherepack_c)
add_subdirectory(spherepack_c)
endif ()

add_subdirectory_if_exists(linal)

if (NOT LINAL_HOME)
        message(FATAL "linal library not found")
endif ()

add_subdirectory(src)
add_subdirectory(utils)

