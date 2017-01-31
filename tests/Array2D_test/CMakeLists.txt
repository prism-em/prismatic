cmake_minimum_required(VERSION 3.5)
project(2D_test)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_CXX_STANDARD 11)
set(Array2D_test_cpu_SOURCE_FILES Array2D_test_construct.cpp)
add_executable(Array2D_test_cpu ${Array2D_test_cpu_SOURCE_FILES})

find_package(CUDA REQUIRED)

set(NVCC_FLAGS_EXTRA "-arch=sm_35 -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})

set(Array2D_test_gpu_SOURCE_FILES Array2D_test_construct.cu)
cuda_include_directories(${CUDA_INCLUDE_DIRS})
cuda_add_executable(Array2D_test_gpu ${Array2D_test_gpu_SOURCE_FILES})

set(Array2D_test_rangefor_cpu_SOURCE_FILES Array2D_test_rangefor.cpp)
cuda_add_executable(Array2D_test_rangefor_cpu ${Array2D_test_rangefor_cpu_SOURCE_FILES})

set(Array2D_test_rangefor_gpu_SOURCE_FILES Array2D_test_rangefor.cu)
cuda_add_executable(Array2D_test_rangefor_gpu ${Array2D_test_rangefor_gpu_SOURCE_FILES})
