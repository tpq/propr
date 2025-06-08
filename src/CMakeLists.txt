add_subdirectory(dispatch)
add_subdirectory(kernels)


set(ALL_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PROPR_SOURCES}
        ${ALL_SOURCES}
)

add_library(propr SHARED ${ALL_SOURCES})
target_link_libraries(propr  ${LIBS} CUDA::cublas)

set_target_properties(propr
        PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${R_LIB}"
            IMPORTED_LOCATION ${RCPP_LIB}
)
target_compile_definitions(propr PUBLIC CONFIG_FILES_DIR="${CMAKE_CONFIG_PREFIX}")
