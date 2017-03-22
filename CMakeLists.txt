cmake_minimum_required(VERSION 3.5)
project(PRISM)
#set(Qt5Widgets_DIR /usr/local/Cellar/qt5/5.8.0_1/lib/cmake/Qt5Widgets/)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O0 -g")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O0")
# PRISM custom options
# Qt options
set(PRISM_ENABLE_GUI 1 CACHE BOOL PRISM_ENABLE_GUI)
set(PRISM_ENABLE_GPU 1 CACHE BOOL PRISM_ENABLE_GPU)
set(PRISM_ENABLE_DOUBLE_PRECISION 0 cache PRISM_ENABLE_DOUBLE_PRECISION)

set(SOURCE_FILES
        src/configure.cpp
        src/getWorkID.cpp
        src/Multislice.cpp
        src/PRISM01.cpp
        src/PRISM02.cpp
        src/PRISM03.cpp
        src/PRISM_entry.cpp
        src/Multislice_entry.cpp
        src/projectedPotential.cpp)

set(GUI_SOURCE_FILES Qt/main.cpp
        Qt/prismmainwindow.cpp)

if (PRISM_ENABLE_GPU)
    set(CUDA_SOURCE_FILES
            src/Multislice.cu
            src/PRISM02.cu
            src/PRISM03.cu
            src/utility.cu)
endif(PRISM_ENABLE_GPU)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

find_package (Threads REQUIRED)
find_package (Boost REQUIRED)
find_package(FFTW REQUIRED)
#find_package(Qt5::Widgets REQUIRED)

set(PRISM_ENABLE_DOUBLE_PRECISION 1)
#set(PRISM_ENABLE_DOUBLE_PRECISION 0)

if (NOT PRISM_ENABLE_DOUBLE_PRECISION)
    message("single precision enabled")
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} fftw3f fftw3f_threads)
else(NOT PRISM_ENABLE_DOUBLE_PRECISION)
    message("double precision enabled")
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} fftw3 fftw3_threads)
endif(NOT PRISM_ENABLE_DOUBLE_PRECISION)

if (PRISM_ENABLE_DOUBLE_PRECISION)
    message("defining double precision macro")
    add_definitions(-DPRISM_ENABLE_DOUBLE_PRECISION)
endif(PRISM_ENABLE_DOUBLE_PRECISION)

#set (PRISM_ENABLE_GPU 0 )
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
endif(PRISM_ENABLE_GUI)
if (PRISM_ENABLE_GUI AND NOT PRISM_ENABLE_GPU)
    	add_executable(prism-gui ${SOURCE_FILES} ${GUI_SOURCE_FILES})
    	target_link_libraries(prism-gui Qt5::Widgets ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
endif (PRISM_ENABLE_GUI AND NOT PRISM_ENABLE_GPU)

if (PRISM_ENABLE_GPU)

    message(GPU enabled)
    set(NVCC_FLAGS_EXTRA " -arch=sm_35 -std=c++11")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})
    cuda_add_library(prism_shared ${SOURCE_FILES} ${CUDA_SOURCE_FILES} SHARED)
    cuda_add_executable(prism src/driver.cpp ${SOURCE_FILES} ${CUDA_SOURCE_FILES})
    cuda_add_cufft_to_target(prism prism_shared)
    include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES} ${CMAKE_SOURCE_DIR})
    target_link_libraries(prism ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
    target_link_libraries(prism_shared ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
    cuda_add_cufft_to_target(prism_shared)
  #  message("FFTW found.   Location " ${FFTW_INCLUDES})
  #  message("Boost found.  Location " ${Boost_Library_DIRS})

if (PRISM_ENABLE_GUI)
    cuda_add_executable(prism-gui  ${SOURCE_FILES} ${CUDA_SOURCE_FILES} ${GUI_SOURCE_FILES} )
    cuda_add_cufft_to_target(prism-gui)

    include_directories(${CMAKE_SOURCE_DIR}/Qt)
    target_link_libraries(prism-gui Qt5::Widgets ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES} prism_shared)
endif(PRISM_ENABLE_GUI)
else()
    add_executable(prism src/driver.cpp ${SOURCE_FILES})
    include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES})
    target_link_libraries(prism ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
    message("FFTW found.   Location " ${FFTW_INCLUDES})
    message("Boost found.  Location " ${Boost_Library_DIRS})
endif (PRISM_ENABLE_GPU)





