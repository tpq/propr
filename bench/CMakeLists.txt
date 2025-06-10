set(BENCHMARK_SOURCES
        lrv.cu
)
foreach(SOURCE IN LISTS BENCHMARK_SOURCES)
    get_filename_component(BENCHMARK_NAME ${SOURCE} NAME_WLE)
    add_executable(${BENCHMARK_NAME} ${SOURCE})
    target_link_libraries(${BENCHMARK_NAME}
            PRIVATE propr
            PRIVATE nvbench::nvbench
    )
    target_include_directories(${BENCHMARK_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/inst/include
    )
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/nvbench)
    get_target_property(PROPR_ARCHITECTURES
            propr CUDA_ARCHITECTURES
    )
    set_target_properties(${BENCHMARK_NAME}
            PROPERTIES
            CUDA_ARCHITECTURES ${PROPR_ARCHITECTURES}
    )
    message(STATUS "Benchmark Added: ${BENCHMARK_NAME}")
endforeach()