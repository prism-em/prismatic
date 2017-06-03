cmake_minimum_required(VERSION 3.5)
project(PRISM)
set(CMAKE_CXX_STANDARD 11)

# set warning level
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w -O3 -fPIC")
endif ( CMAKE_COMPILER_IS_GNUCC )
if ( MSVC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W0 /utf-8 /Ox")
endif ( MSVC )

# PRISM custom options
set(PRISM_ENABLE_GUI 0 CACHE BOOL PRISM_ENABLE_GUI)
set(PRISM_ENABLE_GPU 0 CACHE BOOL PRISM_ENABLE_GPU)
set(PRISM_ENABLE_CLI 1 CACHE BOOL PRISM_ENABLE_GPU)
set(PRISM_ENABLE_DOUBLE_PRECISION 0 CACHE BOOL PRISM_ENABLE_DOUBLE_PRECISION)

#set (CMAKE_BUILD_TYPE DEBUG)
if (PRISM_ENABLE_GPU)
    set(NVCC_FLAGS_EXTRA "-arch=sm_30 -std=c++11 -Xcompiler -fPIC")
endif (PRISM_ENABLE_GPU)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# create variables with various bundles of source files
set(SOURCE_FILES
        src/configure.cpp
        src/WorkDispatcher.cpp
        src/Multislice_calcOutput.cpp
        src/PRISM01_calcPotential.cpp
        src/PRISM02_calcSMatrix.cpp
        src/PRISM03_calcOutput.cpp
        src/PRISM_entry.cpp
        src/Multislice_entry.cpp
	    src/projectedPotential.cpp
        src/parseInput.cpp
        src/utility.cpp
        src/atom.cpp)

if (PRISM_ENABLE_GUI)
set(GUI_SOURCE_FILES
        Qt/main.cpp
        Qt/prismmainwindow.cpp
	Qt/prism_qthreads.cpp
        Qt/prism_progressbar.cpp
#        src/parseInput.cpp
	Qt/prism_resources.qrc)
endif (PRISM_ENABLE_GUI)

if (PRISM_ENABLE_GPU)
    set(CUDA_SOURCE_FILES
            src/Multislice_calcOutput.cu
            src/PRISM02_calcSMatrix.cu
            src/PRISM03_calcOutput.cu
            src/utility.cu)
endif (PRISM_ENABLE_GPU)

# find core packages
find_package (Threads REQUIRED)
find_package (Boost REQUIRED)
find_package(FFTW3 REQUIRED)

message("Boost dir " ${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include
                    ${Boost_INCLUDE_DIRS}
                    ${FFTW_INCLUDE_DIR})
# find CUDA package, if building for GPU
if (PRISM_ENABLE_GPU)
    message("GPU support enabled")
    find_package(CUDA REQUIRED)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${NVCC_FLAGS_EXTRA})
endif(PRISM_ENABLE_GPU)

# find Qt modules, if building GUI
if (PRISM_ENABLE_GUI)
find_package(Qt5Widgets REQUIRED)
message(Qt5Widgets_DIR = ${Qt5Widgets_DIR})
endif (PRISM_ENABLE_GUI)

# add preprocessor directives based on configuration
if (PRISM_ENABLE_DOUBLE_PRECISION)
add_definitions(-DPRISM_ENABLE_DOUBLE_PRECISION)
endif (PRISM_ENABLE_DOUBLE_PRECISION)

if (PRISM_ENABLE_GUI)
    add_definitions(-DPRISM_ENABLE_GUI)
endif (PRISM_ENABLE_GUI)

if (PRISM_ENABLE_GPU)
    add_definitions(-DPRISM_ENABLE_GPU)
endif (PRISM_ENABLE_GPU)

if (PRISM_ENABLE_CLI)
    add_definitions(-DPRISM_ENABLE_CLI)
endif (PRISM_ENABLE_CLI)

# setup FFTW based on desired precision
if (NOT PRISM_ENABLE_DOUBLE_PRECISION)
    message("Single precision enabled")
    if (UNIX)
	    message(UNIX)
    set(FFTW_LIBRARIES ${FFTW_EXECUTABLE_LIBRARIES})
    endif (UNIX)
else(NOT PRISM_ENABLE_DOUBLE_PRECISION)
    message("Double precision enabled")
    add_definitions(-DPRISM_ENABLE_DOUBLE_PRECISION)
	if (FFTW_THREADS_LIBRARY)
    		set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTW_THREADS_LIBRARY})
    	endif(FFTW_THREADS_LIBRARY)
endif(NOT PRISM_ENABLE_DOUBLE_PRECISION)



if (PRISM_ENABLE_CLI)
    # build CLI
    if (PRISM_ENABLE_GPU)

        cuda_add_executable(prism
                src/driver.cpp
                ${SOURCE_FILES}
                ${CUDA_SOURCE_FILES})
        cuda_add_cufft_to_target(prism)
    else(PRISM_ENABLE_GPU)

        add_executable(prism
                        src/driver.cpp
                        ${SOURCE_FILES})
    endif (PRISM_ENABLE_GPU)

    target_link_libraries(prism
        ${CMAKE_THREAD_LIBS_INIT}
        ${Boost_LIBRARY_DIRS}
        ${FFTW_LIBRARIES})
endif (PRISM_ENABLE_CLI)

if (PRISM_ENABLE_GUI)
    # build GUI
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR on)
    include_directories(${CMAKE_SOURCE_DIR}/Qt ${CMAKE_SOURCE_DIR}/include)
    include_directories(${Qt5Widgets_INCLUDE_DIRS})

    if (PRISM_ENABLE_GPU)
        # build the shared CUDA library used by the GUI
	#cuda_add_library(prism_shared
	#                    ${SOURCE_FILES}
	#                    ${CUDA_SOURCE_FILES}
	#                    #${GUI_SOURCE_FILES}
	#                    SHARED)
	#target_link_libraries(prism_shared
	#                        ${CMAKE_THREAD_LIBS_INIT}
	#                        ${Boost_LIBRARY_DIRS}
	#                        ${FFTW_LIBRARIES})
	#cuda_add_cufft_to_target(prism_shared)
#            add_definitions(-DPRISM_BUILDING_GUI)
#        add_definitions(-DHEY)
#
        cuda_add_executable(prism-gui
                                ${SOURCE_FILES}
                                ${GUI_SOURCE_FILES}
                                ${CUDA_SOURCE_FILES}
                                OPTIONS -DPRISM_BUILDING_GUI)
        target_link_libraries(prism-gui
                                Qt5::Widgets
                                ${CMAKE_THREAD_LIBS_INIT}
                                ${Boost_LIBRARY_DIRS}
                                ${FFTW_LIBRARIES})
        cuda_add_cufft_to_target(prism-gui)
    else (PRISM_ENABLE_GPU)
        add_executable(prism-gui
                        ${SOURCE_FILES}
                        ${GUI_SOURCE_FILES})

        target_link_libraries(prism-gui
                                Qt5::Widgets
                                ${CMAKE_THREAD_LIBS_INIT}
                                ${Boost_LIBRARY_DIRS}
                                ${FFTW_LIBRARIES})
    endif (PRISM_ENABLE_GPU)
    target_compile_definitions(prism-gui PRIVATE PRISM_BUILDING_GUI=1)

    install(TARGETS prism-gui RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
endif (PRISM_ENABLE_GUI)

if (PRISM_ENABLE_CLI)
install(TARGETS prism RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
endif(PRISM_ENABLE_CLI)
