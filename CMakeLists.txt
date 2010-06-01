project(sph_ex)
cmake_minimum_required(VERSION 2.8)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CheckedAddSubdir)

find_package(OpenMP)
add_definitions(${OpenMP_C_FLAGS})

add_subdirectory_if_exists(spherepack_c)
add_subdirectory_if_exists(linal)

if (NOT LINAL_HOME)
        message(FATAL "linal library not found")
endif ()

if (NOT SPHEREPACK_INCLUDE)
	message(FATAL "linal spherepack library not found")
endif ()

add_subdirectory(src)
add_subdirectory(utils)

