
project(2D_test)
include_directories(../)
set(CMAKE_CXX_STANDARD 11)
set(TEST_ARRAY_SOURCE_FILES ../Array2D.h Array2D_tester.cpp)
add_executable(2D_test ${TEST_ARRAY_SOURCE_FILES})

