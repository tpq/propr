set(TEST_FILES test_main.cpp)


include(CTest)
add_test(NAME propr-tests COMMAND propr-tests)

enable_testing()
add_executable(propr-tests ${TEST_FILES})
target_include_directories(
        propr-tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/catch2
)
target_link_libraries(propr-tests propr ${LIBS})
target_compile_definitions(propr-tests PUBLIC
                           TEST_FILES_DIR="${CMAKE_TEST_PREFIX}")