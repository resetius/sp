project(examples C)

include_directories(${SPHEREPACK_INCLUDE} .)

if (NOT WIN32)
set(MATH m)
endif ()

set(EXECUTABLE_OUTPUT_PATH "../bin")


add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -DYY_NO_UNISTD_H)
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

flex_target(generator_scanner generator_scanner.l 
	${CMAKE_CURRENT_SOURCE_DIR}/generator_scanner.cpp)
bison_target(generator_parser generator_parser.ypp 
	${CMAKE_CURRENT_SOURCE_DIR}/generator_parser.cpp)
add_flex_bison_dependency(generator_scanner generator_parser)
set(generator_source 
	generator.cpp
	generator.h
	generator_harmonic.cpp
	generator_harmonic.h
	generator_parser.ypp 
	generator_scanner.l 

	${BISON_generator_parser_OUTPUTS}
	${FLEX_generator_scanner_OUTPUTS}
	)
add_library(generator ${generator_source})
add_executable(gen_scheme generator_main.cpp)
target_link_libraries(gen_scheme generator)

add_executable(helm helm.c)
target_link_libraries(helm spherepack ${MATH})

add_executable(helm2 helm.cpp)
target_link_libraries(helm2 spherepack_ext spherepack ${MATH})

add_executable(sp_test_lapl test_lapl.cpp)
target_link_libraries(sp_test_lapl spherepack_ext spherepack ${MATH})

add_executable(sp_test_jac test_jac.cpp)
target_link_libraries(sp_test_jac spherepack_ext spherepack ${MATH})

add_executable(sp_test_chafe test_chafe.cpp)
target_link_libraries(sp_test_chafe spherepack_ext spherepack ${MATH})

add_executable(sp_test_barvortex test_barvortex.cpp)
target_link_libraries(sp_test_barvortex spherepack_ext spherepack ${MATH})

if (OPENMP_FOUND)
	set_target_properties(helm PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(helm2 PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_lapl PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_jac PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_chafe PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_barvortex PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(gen_scheme PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif ()

