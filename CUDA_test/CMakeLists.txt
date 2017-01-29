cmake_minimum_required(VERSION 3.6)
project(CUDA_test)

set(CMAKE_CXX_STANDARD 11)

find_package(CUDA REQUIRED)
if (CUDA_FOUND)
    message("CUDA was found")
endif(CUDA_FOUND)


set(SOURCE_FILES0 CUDA_test1.cu)
add_executable(CUDA_test0 ${SOURCE_FILES0})
set_target_properties(CUDA_test0 PROPERTIES LINKER_LANGUAGE CXX)

set(SOURCE_FILES1 CUDA_test1.cu)
add_executable(CUDA_test1 ${SOURCE_FILES1})
set_target_properties(CUDA_test1 PROPERTIES LINKER_LANGUAGE CXX)

