project(examples C)

include_directories(../spherepack_c)


if (NOT WIN32)
set(MATH m)
endif ()

set(EXECUTABLE_OUTPUT_PATH "../bin")


add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
add_library(spherepack_ext 
	statistics.h 

	grad.cpp 
	grad.h
	lapl.cpp 
	lapl.h
	div.cpp
	div.h
	jac.cpp
	jac.h
	norm.cpp
	norm.h
	vorticity.h
	vorticity.cpp

	chafe.cpp 
	chafe.h
	barvortex.cpp 
	barvortex.h

	utils.cpp
	utils.h
	srtm_rel.cpp
	srtm_rel.h
	)

add_executable(helm helm.c)
target_link_libraries(helm spherepack ${MATH})

add_executable(helm2 helm.cpp)
target_link_libraries(helm2 spherepack_ext spherepack ${MATH})

add_executable(test_lapl test_lapl.cpp)
target_link_libraries(test_lapl spherepack_ext spherepack ${MATH})

add_executable(test_jac test_jac.cpp)
target_link_libraries(test_jac spherepack_ext spherepack ${MATH})

add_executable(test_chafe test_chafe.cpp)
target_link_libraries(test_chafe spherepack_ext spherepack ${MATH})

add_executable(test_barvortex test_barvortex.cpp)
target_link_libraries(test_barvortex spherepack_ext spherepack ${MATH})

