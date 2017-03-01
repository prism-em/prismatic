cmake_minimum_required(VERSION 3.5)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O0")

set(SOURCE_FILES driver.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
add_executable(PRISM ${SOURCE_FILES})
find_package (Threads REQUIRED)
find_package (Boost REQUIRED)
find_package(FFTW REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES})
target_link_libraries(PRISM ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
message("FFTW found.   Location " ${FFTW_INCLUDES})
message("Boost found.  Location " ${FFTW_LIBRARIES})

# PRISM custom options
set(PRISM_ENABLE_GPU 0 CACHE BOOL PRISM_ENABLE_GPU)
if (PRISM_ENABLE_GPU)
    message("GPU support enabled")
endif(PRISM_ENABLE_GPU)

set(PRISM_BUILD_TESTS 0 CACHE BOOL PRISM_BUILD_TESTS)
if (PRISM_BUILD_TESTS)
    add_subdirectory(tests)
    message("Building test programs")
endif(PRISM_BUILD_TESTS)
