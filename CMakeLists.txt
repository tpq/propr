cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0048 NEW)

project(PROPR CXX)

set(PROPR_ENABLE_CUDA "AUTO" CACHE STRING "Enable CUDA backend: AUTO, ON, or OFF")
set_property(CACHE PROPR_ENABLE_CUDA PROPERTY STRINGS AUTO ON OFF)
set(PROPR_HAS_CUDA OFF)

if ("${CMAKE_BUILD_TYPE}" STREQUAL "" OR ${CMAKE_BUILD_TYPE} STREQUAL "NOMODE")
    message("WORKING ON NO MODE")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
    message("WORKING ON RELEASE MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -g")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE} -g")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG" OR ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message("WORKING ON DEBUG MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG}")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
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

if (NOT (PROPR_ENABLE_CUDA STREQUAL "AUTO" OR PROPR_ENABLE_CUDA STREQUAL "ON" OR PROPR_ENABLE_CUDA STREQUAL "OFF"))
    message(FATAL_ERROR "PROPR_ENABLE_CUDA must be one of AUTO, ON, OFF")
endif ()

if (PROPR_ENABLE_CUDA STREQUAL "ON" OR PROPR_ENABLE_CUDA STREQUAL "AUTO")
    check_language(CUDA)
    if (CMAKE_CUDA_COMPILER)
        set(USE_CUDA ON)
        set(PROPR_HAS_CUDA ON)
        message("Build CUDA Support")
        enable_language(CUDA)
        include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --expt-relaxed-constexpr -lineinfo -O3 --use_fast_math --extra-device-vectorization")
        message("CUDA           : FOUND")
        find_package(CUDAToolkit REQUIRED)

        execute_process(
            COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
            COMMAND head -n 1
            COMMAND tr -d .
            OUTPUT_VARIABLE NV_CC
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE NVSIMI_RC
        )

        if (NVSIMI_RC EQUAL 0 AND NOT NV_CC STREQUAL "")
            set(CMAKE_CUDA_ARCHITECTURES ${NV_CC})
            message("Detected Compute Capability : " ${NV_CC})
        else ()
            set(CMAKE_CUDA_ARCHITECTURES 75)
            message("Detected Compute Capability : unavailable, defaulting to 75")
        endif ()
    elseif (PROPR_ENABLE_CUDA STREQUAL "ON")
        message(FATAL_ERROR "CUDA requested (PROPR_ENABLE_CUDA=ON) but CUDA toolkit/compiler not found.")
    else ()
        message("CUDA           : NOT FOUND (building CPU-only backend)")
    endif ()
else ()
    message("CUDA           : DISABLED (PROPR_ENABLE_CUDA=OFF)")
endif ()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr")
endif ()
message("Install Prefix : " ${CMAKE_INSTALL_PREFIX})
include_directories(${CMAKE_SOURCE_DIR}/inst/include)

# set(CUTLASS_NVCC_ARCHS "90")
# add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/cutlass)
# include_directories(${CMAKE_SOURCE_DIR}/thirdparty/cutlass/include/)

option(PROPR_BUILD_TESTS "Option to enable building tests" OFF)

add_subdirectory(src)
message("Install Path   : ${CMAKE_INSTALL_PREFIX}")
message("PROPR Tests    : ${PROPR_BUILD_TESTS}")
message("PROPR CUDA     : ${PROPR_HAS_CUDA}")

install(DIRECTORY inst/include DESTINATION ${CMAKE_INSTALL_PREFIX}/PROPR
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "*.h"
            PATTERN "*.cuh")
