cmake_minimum_required(VERSION 3.6)
project(CUDA_test)

set(CMAKE_CXX_STANDARD 11)

find_package(CUDA REQUIRED)
set(NVCC_FLAGS_EXTRA "-arch=sm_35 -std=c++11")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})

set(SOURCE_FILES0 CUDA_test0.cu)
cuda_add_executable(CUDA_test0 ${SOURCE_FILES0})

set(SOURCE_FILES1 CUDA_test1.cu)
cuda_add_executable(CUDA_test1 ${SOURCE_FILES1})

