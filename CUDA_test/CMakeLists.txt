cmake_minimum_required(VERSION 3.6)
project(CUDA_test C CXX)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES CUDA_test1.cu)
add_executable(CUDA_test ${SOURCE_FILES})
set_target_properties(CUDA_test PROPERTIES LINKER_LANGUAGE CXX)


