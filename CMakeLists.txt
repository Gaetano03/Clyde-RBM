# Minimum version of cmake requried
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# utopiae-orbit project name and languages
project(clyde-rbm LANGUAGES C CXX VERSION 0.0.1)

#if(NOT CMAKE_VERBOSE_MAKEFILE)
set(CMAKE_VERBOSE_MAKEFILE false)
#endif()

# Version of C++ used (C++11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set (CMAKE_BUILD_TYPE Release)
    message(STATUS "Building release version!")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (CMAKE_BUILD_TYPE Debug)
    message (STATUS "WARNING: building debug version!")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -pg -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -pg -fsanitize=undefined")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif ()

#IF (WIN32)
#    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#ELSE ()
#    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
#ENDIF ()

#set(SOPLEX_ROOT_DIR ~/Applications/scipoptsuite-5.0.1/soplex)
#set(SCIP_ROOT_DIR ~/Applications/scipoptsuite-5.0.1/scip)
#find_package(Scip REQUIRED)
#find_package(Threads REQUIRED)

# set variables
set (PROJECT_PATH                                "${CMAKE_CURRENT_SOURCE_DIR}")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}      "${PROJECT_PATH}/cmake/Modules/")
set (INCLUDE_PATH                                "${PROJECT_PATH}/src/")
#set (LIB_PATH                                   "${PROJECT_BINARY_DIR}/lib")
#set (LIB_NAME                                    "clyde-rbm")
#set (LIB_NAME_STATIC                             "clyde-rbm-static")
#set (TEST_PATH                                  "${PROJECT_BINARY_DIR}/test")

# options
option(BUILD_STATIC        "Build static library"     OFF)
option(BUILD_EXAMPLES      "Build examples"           ON)
option(ENABLE_SMARTMATH    "MATH methods"             ON)
option(ENABLE_SMARTUQ      "UQ methods"               ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(GLS REQUIRED)

# Include header directories
get_filename_component(HEADERS "include" REALPATH)
include_directories(${HEADERS} ${GSL_INCLUDE_DIRS})


# Include external libraries
set (MANDATORY_LIBRARIES ${GSL_LIBRARIES})
include(cmake/Dependencies.cmake)




# Include all cpp in src
file(GLOB_RECURSE SOURCES "src/*.cpp")
if (BUILD_STATIC)
    add_library (clyde-rbm STATIC ${SOURCES})
else()
    add_library (clyde-rbm SHARED ${SOURCES})
endif()
set_target_properties(clyde-rbm PROPERTIES LINKER_LANGUAGE CXX)



# Add main as execturable
add_executable(clyde-rbm-main src/main.cpp)
target_link_libraries(clyde-rbm-main clyde-rbm ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(clyde-rbm-main PROPERTIES LINKER_LANGUAGE CXX)
add_dependencies(clyde-rbm-main clyde-rbm)

# If build examples add examples as executables
if ( BUILD_EXAMPLES )
    #set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_PATH})
    add_subdirectory ("examples")
endif()



# Install header files and library.
# Destination is set by CMAKE_INSTALL_PREFIX and defaults to usual locations, unless overridden by
# user.
install (DIRECTORY "${INCLUDE_PATH}/" DESTINATION include/${PROJECT_NAME} FILES_MATCHING PATTERN "*.h")
install (TARGETS clyde-rbm DESTINATION lib/${CMAKE_PROJECT_NAME})
install (TARGETS clyde-rbm DESTINATION bin/${CMAKE_PROJECT_NAME})
