cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(clyde-rbm LANGUAGES C CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_VERBOSE_MAKEFILE FALSE)

if (NOT CMAKE_VERBOSE_MAKEFILE)
    set(CMAKE_VERBOSE_MAKEFILE false)
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS_DEBUG "-g -ggdb -pg -DDEBUG -fsanitize=undefined")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -pg -DDEBUG -fsanitize=undefined")
    set(CMAKE_C_FLAGS_RELEASE "-O4 -DNDEBUG -march=native -fsanitize=undefined")
    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG -march=native -fsanitize=undefined")
endif ()

set(default_build_type "RELEASE")
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "DEBUG")
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "DEBUG" "RELEASE")
endif ()

IF (WIN32)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ELSE ()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF ()

find_package(Threads REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

find_package(GSL REQUIRED)

# Include header directories
get_filename_component(HEADERS "include" REALPATH)
include_directories(${HEADERS} ${GSL_INCLUDE_DIRS})

# Include external libraries
set(MANDATORY_LIBRARIES ${GSL_LIBRARIES})
include(cmake/Dependencies.cmake)

include_directories(${SMART_UQ_INCLUDE_DIR})

# Include all cpp in src
file(GLOB_RECURSE SOURCES "src/*.cpp")
if (BUILD_STATIC)
    add_library(clyde-rbm STATIC ${SOURCES})
else ()
    add_library(clyde-rbm SHARED ${SOURCES})
endif ()
set_target_properties(clyde-rbm PROPERTIES LINKER_LANGUAGE CXX)

if (SMART_UQ_EXTERNAL)
    add_dependencies(clyde-rbm smart-uq-external)
endif ()

# Add main as execturable
add_executable(clyde-rbm-main src/main.cpp)
target_link_libraries(clyde-rbm-main clyde-rbm ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(clyde-rbm-main PROPERTIES LINKER_LANGUAGE CXX)
add_dependencies(clyde-rbm-main clyde-rbm)

# If build examples add examples as executables
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_PATH})
add_subdirectory("examples")
