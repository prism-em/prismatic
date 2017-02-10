cmake_minimum_required(VERSION 3.5)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp include/Array3D.h)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
add_executable(PRISM ${SOURCE_FILES})

find_package(FFTW REQUIRED)
include_directories(FFTW_INCLUDES)
target_link_libraries(PRISM ${FFTW_LIBRARIES})
message("FFTW FOUND! Location " ${FFTW_INCLUDES})

# PRISM custom options
set(PRISM_ENABLE_GPU 1 CACHE BOOL PRISM_ENABLE_GPU)
if (PRISM_ENABLE_GPU)
    message("GPU support enabled")
endif(PRISM_ENABLE_GPU)

set(PRISM_BUILD_TESTS 1 CACHE BOOL PRISM_BUILD_TESTS)
if (PRISM_BUILD_TESTS)
    add_subdirectory(tests)
    message("Building test programs")
endif(PRISM_BUILD_TESTS)