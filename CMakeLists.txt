cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0048 NEW)

project(PROPR CXX)


if ("${CMAKE_BUILD_TYPE}" STREQUAL "" OR ${CMAKE_BUILD_TYPE} STREQUAL "NOMODE")
    message("WORKING ON NO MODE")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
    message("WORKING ON RELEASE MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG" OR ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message("WORKING ON DEBUG MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG}")
else ()
    message(FATAL_ERROR "Unrecognized build type")
endif ()

include(CheckLanguage)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(R REQUIRED)
if (${R_FOUND})
    message("R              : FOUND")
    list(APPEND LIBS "R")
else ()
    message("R              : NOT FOUND")
endif ()

check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    set(USE_CUDA ON)
    message("Build CUDA Support")
    set(PROPERTIES CUDA_ARCHITECTURES 90)
    enable_language(CUDA)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    set(USE_CUDA ON)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}  --extended-lambda")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
    message("CUDA           : FOUND")
    find_package(CUDAToolkit REQUIRED)
    #find_package(Thrust REQUIRED)
else ()
    message("CUDA           : NOT FOUND")
endif ()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr")
endif ()
message("Install Prefix : " ${CMAKE_INSTALL_PREFIX})
include_directories(${CMAKE_SOURCE_DIR}/inst/include)

option(PROPR_BUILD_TESTS "Option to enable building tests" OFF)

add_subdirectory(src)
message("Install Path   : ${CMAKE_INSTALL_PREFIX}")
message("PROPR Tests    : ${PROPR_BUILD_TESTS}")

install(DIRECTORY inst/include DESTINATION ${CMAKE_INSTALL_PREFIX}/PROPR
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "*.h"
            PATTERN "*.cuh")