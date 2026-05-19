set(PROPR_SOURCES "")

add_subdirectory(dispatch)

if (PROPR_BUILD_TESTS)
    find_package(Testthat REQUIRED)
    add_subdirectory(cpp_tests)
endif()

file(GLOB SRC_ROOT_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cu"
)

set(CUDA_PROFILER_REAL "${CMAKE_CURRENT_SOURCE_DIR}/cuda_profiler.cpp")
set(CUDA_PROFILER_STUB "${CMAKE_CURRENT_SOURCE_DIR}/cuda_profiler_stub.cpp")
list(REMOVE_ITEM SRC_ROOT_FILES ${CUDA_PROFILER_REAL} ${CUDA_PROFILER_STUB})

if (PROPR_HAS_CUDA)
    set(CUDA_PROFILER_SOURCE ${CUDA_PROFILER_REAL})
else ()
    set(CUDA_PROFILER_SOURCE ${CUDA_PROFILER_STUB})
endif ()

set(ALL_PROPR_SOURCES ${SRC_ROOT_FILES} ${CUDA_PROFILER_SOURCE} ${PROPR_SOURCES})
add_library(propr SHARED ${ALL_PROPR_SOURCES})
target_link_libraries(propr ${LIBS})
if (PROPR_HAS_CUDA)
    target_link_libraries(propr CUDA::cublas CUDA::cudart)
endif ()
# target_include_directories(propr
#         PRIVATE
#         ${CMAKE_SOURCE_DIR}/thirdparty/cutlass/include
# )

if (PROPR_BUILD_TESTS)
    target_link_libraries(propr Testthat::Testthat)
endif()

if (PROPR_HAS_CUDA AND CMAKE_CUDA_ARCHITECTURES)
    set_property(TARGET propr PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif ()
set_target_properties(propr
        PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${R_LIB}"
            IMPORTED_LOCATION ${RCPP_LIB}
)
target_compile_definitions(propr PUBLIC CONFIG_FILES_DIR="${CMAKE_CONFIG_PREFIX}")
if (PROPR_HAS_CUDA)
    target_compile_definitions(propr PUBLIC PROPR_HAS_CUDA=1)
else ()
    target_compile_definitions(propr PUBLIC PROPR_HAS_CUDA=0)
endif ()
