## CppUTest_FOUND        = CppUTest is installed
## CppUTest_INCLUDE_DIR = CppUTest include directory
## CppUTest_LIBRARIES    = Link options to compile with CppUTest
## CppUTest_LIBRARY      = The CppUTest library
#
## You can use the CppUTest_HOME environment variable
#
## References:
## [1] https://github.com/Adnn/aunteater/tree/master/cmake
## [2] http://noqisofon.hatenablog.com/entry/20120730/1343631914

set(CppUTest_HOME $ENV{CppUTest_HOME} CACHE PATH "Path to CppUTest" )

##
## Look for headers and libs
##
find_path(CppUTest_INCLUDE_DIRS
    NAMES CppUTest/TestHarness.h
    PATHS
        $ENV{CppUTest_HOME}/include
        /usr/local/include
        /usr/include
  )

find_library(CppUTest_LIBRARIES
    NAMES libCppUTest CppUTest
    PATHS
        ${CppUTest_HOME}
        ${CppUTest_HOME}/lib
        ${CppUTest_HOME}/lib/release
        ${CppUTest_HOME}/bin
        /usr/local/lib
        /usr/lib
        /usr/lib64
  )


##
## handle the QUIETLY and REQUIRED arguments and set CppUTest_FOUND to TRUE if
## all listed variables are TRUE
##
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CppUTest DEFAULT_MSG CppUTest_INCLUDE_DIRS CppUTest_LIBRARIES)

MARK_AS_ADVANCED( CppUTest_INCLUDE_DIRS CppUTest_LIBRARIES )
