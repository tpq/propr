set(TEST_SOURCES test_main.cu test_lrv.cu)
add_executable(test_main ${TEST_SOURCES})
target_link_libraries(test_main  CUDA::cublas CUDA::cudart propr Catch2::Catch2)
set_property(TARGET test_main PROPERTY CUDA_ARCHITECTURES 90)
add_test(NAME test_main COMMAND test_main)
target_include_directories(test_main PRIVATE
        ${CMAKE_SOURCE_DIR}/inst/include
        ${CMAKE_CURRENT_SOURCE_DIR}
)
set_target_properties(test_main
        PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE}"
            INTERFACE_LINK_LIBRARIES "${R_LIB}"
            IMPORTED_LOCATION ${RCPP_LIB}
)

#set (CMAKE_CXX_FLAGS "-stdlib=libstdc++ ${CMAKE_CXX_FLAGS}")
#LDFLAGS=-stdlib=libstdc++
#file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cu")
#foreach(test_src IN LISTS TEST_SOURCES)
#    get_filename_component(test_name ${test_src} NAME_WE)
#    add_executable(${test_name} ${test_src})
#    target_link_libraries(${test_name}
#            PRIVATE
#            propr
#            Catch2::Catch2
#    )
#    target_include_directories(${test_name} PRIVATE
#            ${CMAKE_SOURCE_DIR}/inst/include
#            ${CMAKE_CURRENT_SOURCE_DIR}
#    )
#    add_test(NAME ${test_name} COMMAND ${test_name})
#endforeach()
