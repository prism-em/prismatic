cmake_minimum_required(VERSION 3.6)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)
add_subdirectory(Array2D_test)
add_subdirectory(CUDA_test)

set(SOURCE_FILES main.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
add_executable(PRISM ${SOURCE_FILES})

find_package(FFTW REQUIRED)
include_directories(FFTW_INCLUDES)
target_link_libraries(PRISM ${FFTW_LIBRARIES})
message("FFTW FOUND! Location " ${FFTW_INCLUDES})

