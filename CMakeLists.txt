CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(eertools CXX)

SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-deprecated -ftree-vectorize")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_BUILD_TYPE Release)

#if(${CMAKE_COMPILER_IS_GNUCXX})
	add_definitions(-std=c++11)
#endif()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin/")

SET(serial "1.0.0")
SET(soserial "1")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
ADD_LIBRARY(eertools SHARED eertools.cpp)
SET_TARGET_PROPERTIES(eertools PROPERTIES VERSION ${serial} SOVERSION ${soserial})

# Check for CppUTest and bail out if they don't have it
if(DEFINED ENV{CPPUTEST_HOME})
	enable_testing()
	message("Using CppUTest found in $ENV{CPPUTEST_HOME}")

	add_executable(eertools_cpputest eertools_cpputest.cpp)
	target_link_libraries(eertools_cpputest eertools)
	target_link_libraries(eertools_cpputest CppUTest)
	add_test(eertools_cpputest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/eertools_cpputest)
else()
	message("CPPUTEST_HOME is not set; You must tell CMake where to find CppUTest")
	return()
endif()

INSTALL(TARGETS eertools LIBRARY DESTINATION lib)
INSTALL(FILES eertools.h DESTINATION include)
