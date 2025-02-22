enable_testing()
include(GoogleTest)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(test_solver test.cpp)
target_link_libraries(test_solver gtest gtest_main ${PROJECT_NAME})
set_target_properties(test_solver PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})



gtest_discover_tests(test_solver)
add_test(NAME TridiagonalSolverTest COMMAND test_solver)


