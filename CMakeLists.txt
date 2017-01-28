cmake_minimum_required(VERSION 3.6)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
add_executable(PRISM ${SOURCE_FILES})

set(TEST_ARRAY_SOURCE_FILES Array2D.h Array2D_tester.cpp)
add_executable(Array2D_tester ${TEST_ARRAY_SOURCE_FILES})

find_package(FFTW REQUIRED)
include_directories(FFTW_INCLUDES)
target_link_libraries(PRISM ${FFTW_LIBRARIES})
message("FFTW FOUND! Location " ${FFTW_INCLUDES})