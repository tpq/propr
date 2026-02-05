add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cuda)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cpu)
set(PROPR_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/device_selector.cpp
        ${PROPR_SOURCES}
        PARENT_SCOPE)