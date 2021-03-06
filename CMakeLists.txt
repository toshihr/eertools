CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
	if(POLICY CMP0042)
		cmake_policy(SET CMP0042 NEW)
	endif()
endif(COMMAND cmake_policy)

PROJECT(EERTOOLS CXX)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# === VERSION ===
SET(VERSION "1.0.1")
SET(SOVERSION "1")

# === COMPILE OPTIONS ===
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-deprecated -ftree-vectorize")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_BUILD_TYPE Release)

# === C++11 CHECK ===
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# === CPPUTEST ===
find_package(CppUTest)

ADD_SUBDIRECTORY(src)
if(CPPUTEST_FOUND)
	ADD_SUBDIRECTORY(test)
	ENABLE_TESTING()
	ADD_TEST(NAME MyTest COMMAND Test)
endif()


# References
# [1] http://stackoverflow.com/questions/14446495/cmake-project-structure-with-unit-tests
