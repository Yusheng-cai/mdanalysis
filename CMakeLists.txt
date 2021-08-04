cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(mdanalysis)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -Wall")

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

enable_testing()

add_subdirectory(xdr)
add_subdirectory(tools)
add_subdirectory(OrderParameters)
add_subdirectory(parallel)
add_subdirectory(test)

add_executable(mdanalysis main.cpp)

target_include_directories(mdanalysis PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(mdanalysis PUBLIC OrderParameters xdr tools parallel)

install(TARGETS mdanalysis DESTINATION bin)
    