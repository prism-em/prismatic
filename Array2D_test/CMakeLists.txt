cmake_minimum_required(VERSION 3.5)
project(2D_test)
include_directories(../)

set(CMAKE_CXX_STANDARD 11)
set(TEST_ARRAY_SOURCE_FILES ../Array2D.h Array2D_tester.cpp)
add_executable(Array2D_test_cpu ${TEST_ARRAY_SOURCE_FILES})

find_package(CUDA REQUIRED)

set(NVCC_FLAGS_EXTRA "-arch=sm_35 -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})

set(TEST_ARRAY_SOURCE_FILES ../Array2D.h Array2D_tester.cu)
cuda_include_directories(${CUDA_INCLUDE_DIRS})
cuda_add_executable(Array2D_test_gpu ${TEST_ARRAY_SOURCE_FILES})

