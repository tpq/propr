cmake_minimum_required(VERSION 3.15)

set(CMAKE_MODULE_PATH /home/zalbanob/propr/cmake)
set(RCPP_INCLUDE_DIRS /home/zalbanob/R/x86_64-pc-linux-gnu-library/4.5/Rcpp/include)
set(RCPP_LIB_PATH /home/zalbanob/R/x86_64-pc-linux-gnu-library/4.5/Rcpp/libs)

set(CMAKE_CUDA_ARCHITECTURES 80 90)
set(CUDACXX "/usr/local/cuda/bin/nvcc" CACHE FILEPATH "Path to nvcc compiler")
project(propr LANGUAGES CXX CUDA)

# --- Configuration Options ---

# Set C++ standard globally for the project.
# Modern CMake prefers target_compile_features per target for better propagation.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Generally good practice
# --- Find External Dependencies ---

# Find R. This provides R_INCLUDE_DIRS, R_LIBRARIES, R_BINARY etc.
# R_BINARY is crucial for querying R for Rcpp flags.
find_package(R REQUIRED)

# Rcpp does not provide a standard CMake package config in the usual sense.
# We simulate R's Makevars by querying the R executable for Rcpp's CXXFLAGS and LDFLAGS.
# This ensures compatibility with how R packages normally build.

# Get Rcpp CXX flags (primarily for include paths)
execute_process(
        COMMAND "${R_BINARY}" -e "cat(Rcpp:::CxxFlags())"
        OUTPUT_VARIABLE RCPP_CXX_FLAGS_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Parse the -I paths from the output.
# Example output: "-I/path/to/R/library/Rcpp/include"
string(REPLACE "-I" "" RCPP_INCLUDE_DIRS_RAW "${RCPP_CXX_FLAGS_OUTPUT}")
string(REPLACE " " ";" RCPP_INCLUDE_DIRS "${RCPP_INCLUDE_DIRS_RAW}")

# Get Rcpp Linker flags (primarily for library paths and names)
execute_process(
        COMMAND "${R_BINARY}" -e "cat(Rcpp:::LdFlags())"
        OUTPUT_VARIABLE RCPP_LDFLAGS_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Parse the -L (library search path) and -l (library name) parts from RCPP_LDFLAGS_OUTPUT
set(RCPP_LINK_DIRS "")
set(RCPP_LIBRARIES_LIST "")
string(REGEX MATCHALL "-L[^ ]+" RCPP_LIB_DIRS_MATCHES "${RCPP_LDFLAGS_OUTPUT}")
foreach(MATCH ${RCPP_LIB_DIRS_MATCHES})
    string(REPLACE "-L" "" LIB_DIR ${MATCH})
    list(APPEND RCPP_LINK_DIRS ${LIB_DIR})
endforeach()

string(REGEX MATCHALL "-l[^ ]+" RCPP_LIB_NAMES_MATCHES "${RCPP_LDFLAGS_OUTPUT}")
foreach(MATCH ${RCPP_LIB_NAMES_MATCHES})
    string(REPLACE "-l" "" LIB_NAME ${MATCH})
    list(APPEND RCPP_LIBRARIES_LIST ${LIB_NAME})
endforeach()

# Find CUDA Toolkit. This provides CUDA_TOOLKIT_FOUND, CUDAToolkit_INCLUDE_DIRS, CUDAToolkit_LIBRARIES.
# Ensure the version matches what you expect.
find_package(CUDAToolkit 11.0 REQUIRED)


set(PROPR_LIB_SOURCES
        src/core/kernels/cuda/lrm.cu               # CUDA source
        src/core/kernels/cuda/lrv.cu               # CUDA source

        src/core/dispatch/cpu/backend.cpp
        src/core/dispatch/cpu/comparison.cpp
        src/core/dispatch/cpu/ctzRcpp.cpp
        src/core/dispatch/cpu/graflex.cpp
        src/core/dispatch/cpu/lr2propr.cpp
        src/core/dispatch/cpu/lrm.cpp
        src/core/dispatch/cpu/lrv.cpp
        src/core/dispatch/cpu/omega.cpp
        src/core/kernels/cuda/omega.cu-------old             # CUDA source
        src/core/RcppExports.cpp
        src/include/context.h
        src/include/interface/device_selector.hpp
        src/core/dispatch/device_selector.cpp
        src/core/interface/lrv.cpp
        src/core/interface/lrm.cpp
        src/core/dispatch/cpu/lrm.cpp
        src/core/dispatch/cpu/backend.cpp
        src/include/kernels/cpu/dispatch/backend.hpp
        src/core/interface/backend.cpp
        src/core/interface/comparison.cpp
        src/core/interface/graflex.cpp
        src/core/interface/lr2propr.cpp
        src/core/interface/omega.cpp
)

set(PROPR_BENCHMARK_SOURCES
        bench/src/main_benchmark.cu # Benchmark source
)

set(PROPR_TEST_SOURCES
        tests/cpp_test/test_lrm.cpp
        tests/cpp_test/test_lrv.cpp
        tests/cpp_test/test_main.cpp
        tests/cpp_test/test_omega.cpp
)


add_library(propr SHARED ${PROPR_LIB_SOURCES})

target_compile_features(propr PUBLIC cxx_std_14)
set_target_properties(propr PROPERTIES
        CUDA_STANDARD 14
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

# Specify include directories for the library.
target_include_directories(propr PUBLIC
        # R includes (from find_package(R))
        ${R_INCLUDE_DIRS}
        # Rcpp includes (parsed from Rcpp:::CxxFlags())
        ${RCPP_INCLUDE_DIRS}
        # CUDA Toolkit includes (from find_package(CUDAToolkit))
        ${CUDAToolkit_INCLUDE_DIRS}
        # Local project headers (e.g., src/include/interface, src/include/kernels etc.)
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include
        # Prereq headers (e.g., Catch2 for tests, if they are part of the main library)
        # The original CMakeLists.txt had `include_directories(prereq)`.
        # Assuming `prereq` includes `prereq/catch2` and other necessary headers.
        ${CMAKE_CURRENT_SOURCE_DIR}/prereq
)

# Link libraries for the main shared library.
target_link_libraries(propr PUBLIC
        # R runtime libraries (from find_package(R))
        ${R_LIBRARIES}
        # Rcpp libraries (parsed from Rcpp:::LdFlags())
        ${RCPP_LIBRARIES_LIST}
        # CUDA runtime libraries (from find_package(CUDAToolkit))
        ${CUDAToolkit_LIBRARIES}
)

# Add Rcpp library search paths (from -L flags)
target_link_directories(propr PUBLIC ${RCPP_LINK_DIRS})

# R package specific output naming conventions:
# R expects the shared library to be named 'propr.so' (on Linux/macOS)
# or 'propr.dll' (on Windows). CMake typically adds 'lib' prefix.
# We override this for the R package.
set_target_properties(propr PROPERTIES
        OUTPUT_NAME "propr"
        # CMake will automatically use .dll for Windows and .dylib for macOS when building.
        # We specify .so as a common default, but the CMAKE_SHARED_LIBRARY_SUFFIX variable
        # will ultimately determine the actual suffix based on the platform.
        SUFFIX ".so"
)

# 2. Benchmark Executable
# This is a standalone executable for performance testing.
add_executable(propr_benchmark ${PROPR_BENCHMARK_SOURCES})

# Link against the main library (propr) so it can use its functions.
target_link_libraries(propr_benchmark PRIVATE propr)

# Ensure benchmark can find necessary headers. It should implicitly get public includes from 'propr'.
# Add any private includes specific to the benchmark source.
target_include_directories(propr_benchmark PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include # If benchmark uses internal headers not publicly exposed by 'propr'
)

# Set compile features for benchmark.
target_compile_features(propr_benchmark PUBLIC cxx_std_14)
set_target_properties(propr_benchmark PROPERTIES
        CUDA_STANDARD 14
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

# 3. Test Executable
# This is a standalone executable for running unit tests (e.g., with Catch2).
add_executable(propr_tests ${PROPR_TEST_SOURCES})

# Link against the main library (propr) to test its functions.
target_link_libraries(propr_tests PRIVATE propr)

# Add include directories for Catch2 (assuming it's a header-only library in prereq/catch2)
target_include_directories(propr_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/prereq/catch2
        ${CMAKE_CURRENT_SOURCE_DIR}/src/include # If tests use internal headers not publicly exposed by 'propr'
)

# Set compile features for tests.
target_compile_features(propr_tests PUBLIC cxx_std_14)
set_target_properties(propr_tests PROPERTIES
        CUDA_STANDARD 14
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

# Enable CTest for running tests
enable_testing()
add_test(NAME run_propr_cpp_tests COMMAND propr_tests)
