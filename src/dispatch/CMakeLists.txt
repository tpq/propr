add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/runtime)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cpu)

if (PROPR_HAS_CUDA)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cuda)
endif()

set(PROPR_SOURCES
        ${PROPR_SOURCES}
        PARENT_SCOPE)
