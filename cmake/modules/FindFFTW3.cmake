## https://github.com/thierry-sousbie/dice/edit/master/modules/FindFFTW3.cmake
# - Find FFTW
# Find the native FFTW includes and library
# This module defines
# FFTW_INCLUDE_DIR, where to find fftw3.h, etc.
# FFTW_LIBRARIES, the libraries needed to use FFTW.
# FFTW_FOUND, If false, do not try to use FFTW.
# also defined, but not for general use are
# FFTW_LIBRARY, where to find the FFTW library.

#UNSET(FFTW_INCLUDE_DIR CACHE)
#UNSET(FFTW_LIBRARY CACHE)
#UNSET(FFTW_THREADS_NAMES CACHE)
#UNSET(FFTW_MPI_LIBRARY CACHE)

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
        PATHS ${FFTW3_DIR}
        PATH_SUFFIXES include
        DOC "Directory where the FFTW3 header files are located"
        NO_DEFAULT_PATH
        )

SET(FFTW_NAMES ${FFTW_NAMES} fftw3 fftw3-3)
find_library(FFTW_LIBRARY
        NAMES ${FFTW_NAMES}
        PATHS ${FFTW3_DIR}
        PATH_SUFFIXES lib
        DOC "Directory where the FFTW3 library is located"
        NO_DEFAULT_PATH
        )

SET(FFTW_THREADS_NAMES ${FFTW_THREADS_NAMES} fftw3_threads fftw3-3_threads)
find_library(FFTW_THREADS_LIBRARY
        NAMES ${FFTW_THREADS_NAMES}
        PATHS ${FFTW3_DIR}
        PATH_SUFFIXES lib
        DOC "Directory where the FFTW3-threads library is located"
        NO_DEFAULT_PATH
        )

SET(FFTW_MPI_NAMES ${FFTW_MPI_NAMES} fftw3_mpi fftw3-3_mpi)
find_library(FFTW_MPI_LIBRARY
        NAMES ${FFTW_MPI_NAMES}
        PATHS ${FFTW3_DIR}
        PATH_SUFFIXES lib
        DOC "Directory where the FFTW3-MPI library is located"
        NO_DEFAULT_PATH
        )

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
        ${FFTW3_DIR} ${FFTW3_DIR}/include ${FFTW3_DIR}/lib )# /usr/local/include /usr/include /opt/local/lib )

FIND_LIBRARY(FFTW_LIBRARY
        NAMES ${FFTW_NAMES}
        PATHS ${FFTW3_DIR} ${FFTW3_DIR}/include ${FFTW3_DIR}/lib )#/usr/lib /usr/local/lib /opt/locala/lib )

# Find threads part of FFTW
FIND_LIBRARY(FFTW_THREADS_LIBRARY
        NAMES ${FFTW_THREADS_NAMES}
        PATHS ${FFTW3_DIR} ${FFTW3_DIR}/include ${FFTW3_DIR}/lib )#/usr/lib /usr/local/lib /opt/local/lib )

# Find MPI part of FFTW
FIND_LIBRARY(FFTW_MPI_LIBRARY
        NAMES ${FFTW_MPI_NAMES}
        PATHS ${FFTW3_DIR} ${FFTW3_DIR}/include ${FFTW3_DIR}/lib )#/usr/lib /usr/local/lib /opt/local/lib )

IF (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_THREADS_LIBRARIES ${FFTW_THREADS_LIBRARY})
    SET(FFTW_THREADS_FOUND "YES")
ELSE (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_THREADS_FOUND "NO")
ENDIF (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)


IF (FFTW_THREADS_FOUND)
    IF (NOT FFTW_THREADS_FIND_QUIETLY)
        MESSAGE(STATUS "Found FFTW threads: ${FFTW_THREADS_LIBRARIES}")
    ENDIF (NOT FFTW_THREADS_FIND_QUIETLY)
ELSE (FFTW_THREADS_FOUND)
    IF (FFTW_THREADS_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find FFTW threads library")
    ENDIF (FFTW_THREADS_FIND_REQUIRED)
ENDIF (FFTW_THREADS_FOUND)


IF (FFTW_MPI_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_MPI_LIBRARIES ${FFTW_MPI_LIBRARY})
    SET(FFTW_MPI_FOUND "YES")
ELSE (FFTW_MPI_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_MPI_FOUND "NO")
ENDIF (FFTW_MPI_LIBRARY AND FFTW_INCLUDE_DIR)


IF (FFTW_MPI_FOUND)
    IF (NOT FFTW_MPI_FIND_QUIETLY)
        MESSAGE(STATUS "Found FFTW MPI: ${FFTW_MPI_LIBRARIES}")
    ENDIF (NOT FFTW_MPI_FIND_QUIETLY)
ELSE (FFTW_MPI_FOUND)
    IF (FFTW_MPI_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find FFTW MPI library")
    ENDIF (FFTW_MPI_FIND_REQUIRED)
ENDIF (FFTW_MPI_FOUND)


IF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_LIBRARIES ${FFTW_LIBRARY})
    SET(FFTW_FOUND "YES")
ELSE (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_FOUND "NO")
ENDIF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)


IF (FFTW_FOUND)
    IF (NOT FFTW_FIND_QUIETLY)
        MESSAGE(STATUS "Found FFTW: ${FFTW_LIBRARIES}")
    ENDIF (NOT FFTW_FIND_QUIETLY)
ELSE (FFTW_FOUND)
    IF (FFTW_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find FFTW library")
    ENDIF (FFTW_FIND_REQUIRED)
ENDIF (FFTW_FOUND)

SET (ON_UNIX ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR
        ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
IF (${ON_UNIX})
    SET (FFTW_EXECUTABLE_LIBRARIES fftw3f fftw3f_threads)
ENDIF (${ON_UNIX})
