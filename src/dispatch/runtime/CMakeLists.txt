set(RUNTIME_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/resolve_backend.cpp
)

if (NOT PROPR_HAS_CUDA)
    list(APPEND RUNTIME_SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/stubs/cuda_stubs.cpp
    )
endif()

set(PROPR_SOURCES
        ${RUNTIME_SOURCES}
        ${PROPR_SOURCES}
        PARENT_SCOPE
)
