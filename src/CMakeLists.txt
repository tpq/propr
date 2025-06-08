set(PROPR_SOURCES "")

add_subdirectory(dispatch)
add_subdirectory(kernels)

file(GLOB SRC_ROOT_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cu"
)

set(ALL_PROPR_SOURCES ${SRC_ROOT_FILES} ${PROPR_SOURCES})
add_library(propr SHARED ${ALL_PROPR_SOURCES})
target_link_libraries(propr  ${LIBS} CUDA::cublas CUDA::cudart)

set_target_properties(propr
        PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${R_LIB}"
            IMPORTED_LOCATION ${RCPP_LIB}
)
target_compile_definitions(propr PUBLIC CONFIG_FILES_DIR="${CMAKE_CONFIG_PREFIX}")