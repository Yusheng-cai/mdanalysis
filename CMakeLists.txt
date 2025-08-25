cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(mdanalysis)

#find_package(Torch REQUIRED)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require C++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -ffast-math")

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Requires boost
find_package(Boost REQUIRED)
if (Boost_FOUND)
    message("Boost is found")
endif()

# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

# Requires FFTW3 for now 
find_package(FFTW3)
if (FFTW3_FOUND)
  message(STATUS "FFTW3 library: ${FFTW3_LIBRARIES}")
endif()


enable_testing()

add_subdirectory(xdr)
add_subdirectory(tools)
add_subdirectory(OrderParameters)
add_subdirectory(parallel)
add_subdirectory(Actions)
add_subdirectory(test)

add_executable(mdanalysis main.cpp)
add_executable(mda_cmd mda_cmd.cpp)

target_include_directories(mdanalysis PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(mda_cmd PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(mdanalysis PUBLIC mdaOrderParameters mdaxdr mdatools mdaparallel)
target_link_libraries(mda_cmd PUBLIC mdaActions mdaOrderParameters mdaxdr mdatools mdaparallel)

install(TARGETS mdanalysis DESTINATION bin)
install(TARGETS mda_cmd DESTINATION bin)