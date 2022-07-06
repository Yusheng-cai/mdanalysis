# A cmake script that searches for FFTW3 library 
# If we were to include this directory (cmake), when we do 
# FIND_PACKAGE(foo), it will search for the file Findfoo.cmake


INCLUDE(FindPackageHandleStandardArgs)

# Checks the environmental variable
if(FFTW3_DIR)
	# Use the hint provided
	list(APPEND CMAKE_PREFIX_PATH ${FFTW3_DIR})
endif() 

set(CMAKE_FIND_LIBRARY_SUFFIXES .a)

find_path(
    FFTW3_INCLUDE_DIR 
        fftw3.h
    HINTS
        ${FFTW3_DIR}
)

find_library(
    FFTW3_LIBRARY
	NAMES fftw3f
)


find_package_handle_standard_args( FFTW3 DEFAULT_MSG
  FFTW3_INCLUDE_DIR
  FFTW3_LIBRARY
)

if (FFTW3_FOUND)
    SET( FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR} )
    SET( FFTW3_LIBRARIES ${FFTW3_LIBRARY} )
    message("Found FFTW3. FFTW3_LIBRARY: ${FFTW3_LIBRARIES}. FFTW3_INCLUDE_DIRS: ${FFTW3_INCLUDE_DIRS}.")
endif()
