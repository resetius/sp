project(examples C)

include_directories(../spherepack_c)

set(EXECUTABLE_OUTPUT_PATH "../bin")

add_library(spherepack_ext lapl.cpp)

add_executable(helm helm.c)
target_link_libraries(helm spherepack m)

add_executable(helm2 helm.cpp)
target_link_libraries(helm2 spherepack_ext spherepack m)

add_executable(test_lapl test_lapl.cpp)
target_link_libraries(test_lapl spherepack_ext spherepack m)

