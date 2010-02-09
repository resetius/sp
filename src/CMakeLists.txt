project(examples C)

include_directories(../spherepack_c)
add_executable(helm helm.c)
target_link_libraries(helm spherepack m)

