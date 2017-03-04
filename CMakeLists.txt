cmake_minimum_required(VERSION 3.5)
project(PRISM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")

set(SOURCE_FILES driver.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
add_executable(prism ${SOURCE_FILES})
find_package (Threads REQUIRED)
find_package (Boost REQUIRED)
find_package(FFTW REQUIRED)
include_directories(${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDES})
target_link_libraries(prism ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES})
message("FFTW found.   Location " ${FFTW_INCLUDES})
message("Boost found.  Location " ${Boost_Library_DIRS})

# PRISM custom options

# Qt options
set(PRISM_ENABLE_GUI 1 CACHE BOOL PRISM_ENABLE_GUI)
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

set(PRISM_ENABLE_GPU 0 CACHE BOOL PRISM_ENABLE_GPU)
if (PRISM_ENABLE_GPU)
    message("GPU support enabled")
endif(PRISM_ENABLE_GPU)

