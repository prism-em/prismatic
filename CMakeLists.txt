cmake_minimum_required(VERSION 3.5)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")
# PRISM custom options
# Qt options
set(PRISM_ENABLE_GUI 0 CACHE BOOL PRISM_ENABLE_GUI)
set(PRISM_ENABLE_GPU 0 CACHE BOOL PRISM_ENABLE_GPU)
set(PRISM_ENABLE_SINGLE_PRECISION 1 cache PRISM_ENABLE_SINGLE_PRECISION)

if (PRISM_ENABLE_DOUBLE_PRECISION)
    add_definitions(-DPRISM_ENABLE_DOUBLE_PRECISION)
endif(PRISM_ENABLE_DOUBLE_PRECISION)

set(PRISM_ENABLE_GPU 0)
if (PRISM_ENABLE_GPU)
    find_package(CUDA REQUIRED)
    add_definitions(-DPRISM_ENABLE_GPU)
    message("GPU support enabled")
endif(PRISM_ENABLE_GPU)

#set(Qt5Widgets_DIR)
set(CMAKE_INCLUDE_CURRENT_DIR on)
if (PRISM_ENABLE_GUI)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    find_package(Qt5Widgets REQUIRED)
    message(Qt5Widgets_DIR = ${Qt5Widgets_DIR})
    include_directories(${CMAKE_SOURCE_DIR}/Qt)
    add_executable(prism-gui Qt/main.cpp Qt/prismmainwindow.cpp)
    target_link_libraries(prism-gui Qt5::Widgets ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
endif (PRISM_ENABLE_GUI)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

find_package (Threads REQUIRED)
find_package (Boost REQUIRED)
find_package(FFTW REQUIRED)

if (PRISM_ENABLE_SINGLE_PRECISION)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} fftw3f)
endif(PRISM_ENABLE_SINGLE_PRECISION)


set(SOURCE_FILES driver.cpp configure.cpp)
if (PRISM_ENABLE_GPU)
    message(GPU enabled)
    set(NVCC_FLAGS_EXTRA "-arch=sm_35 -std=c++11")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})
    set(SOURCE_FILES ${SOURCE_FILES} include/Multislice.cu)
    cuda_add_executable(prism ${SOURCE_FILES})
    cuda_add_cufft_to_target(prism)
    include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES} ${CMAKE_SOURCE_DIR})
    target_link_libraries(prism ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
    message("FFTW found.   Location " ${FFTW_INCLUDES})
    message("Boost found.  Location " ${Boost_Library_DIRS})
else()
    add_executable(prism ${SOURCE_FILES})
    include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES})
    target_link_libraries(prism ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
    message("FFTW found.   Location " ${FFTW_INCLUDES})
    message("Boost found.  Location " ${Boost_Library_DIRS})
endif (PRISM_ENABLE_GPU)





