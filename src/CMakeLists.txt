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

set(ALL_PROPR_SOURCES ${SRC_ROOT_FILES} ${PROPR_SOURCES})
add_library(propr SHARED ${ALL_PROPR_SOURCES})
target_link_libraries(propr  ${LIBS} CUDA::cublas CUDA::cudart)
target_include_directories(propr
        PRIVATE
        ${CMAKE_SOURCE_DIR}/thirdparty/cutlass/include
)

if (PROPR_BUILD_TESTS)
    target_link_libraries(propr Testthat::Testthat)
endif()

set_property(TARGET propr PROPERTY CUDA_ARCHITECTURES ${NV_CC})
set_target_properties(propr
        PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${R_LIB}"
            IMPORTED_LOCATION ${RCPP_LIB}
)
target_compile_definitions(propr PUBLIC CONFIG_FILES_DIR="${CMAKE_CONFIG_PREFIX}")