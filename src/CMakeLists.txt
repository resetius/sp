project(sp_src C)

include_directories(${SPHEREPACK_INCLUDE} ${LINAL_INCLUDE} .)

if (NOT WIN32)
set(MATH m)
endif ()

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -DYY_NO_UNISTD_H ${LINAL_DEFINES})
add_library(sp 
	statistics.h 

	operator.cpp
	operator.h
	norm.cpp
	norm.h
	grad.cpp 
	grad.h
	lapl.cpp 
	lapl.h
	div.cpp
	div.h
	jac.cpp
	jac.h
	vorticity.h
	vorticity.cpp

	chafe.cpp 
	chafe.h
	barvortex.cpp 
	barvortex.h
	baroclin.cpp
	baroclin.h
	)

add_executable(sp_test_helm helm.c)
target_link_libraries(sp_test_helm spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_helm2 helm.cpp)
target_link_libraries(sp_test_helm2 sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_lapl test_lapl.cpp)
target_link_libraries(sp_test_lapl sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_jac test_jac.cpp)
target_link_libraries(sp_test_jac sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_chafe test_chafe.cpp)
target_link_libraries(sp_test_chafe sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_barvortex test_barvortex.cpp)
target_link_libraries(sp_test_barvortex sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_forcing test_forcing.cpp)
target_link_libraries(sp_test_forcing sp spherepack ${MATH} ${LINAL_LIBS})

add_executable(sp_test_baroclin test_baroclin.cpp)
target_link_libraries(sp_test_baroclin sp spherepack ${MATH} ${LINAL_LIBS})

if (OPENMP_FOUND)
	set_target_properties(sp_test_helm PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_helm2 PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_lapl PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_jac PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_chafe PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_barvortex PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_forcing PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(sp_test_baroclin PROPERTIES LINK_FLAGS ${OpenMP_C_FLAGS})
endif ()

add_test(sp_test_barvortex ${EXECUTABLE_OUTPUT_PATH}/sp_test_barvortex test)

