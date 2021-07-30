add_library(parallel SHARED)

target_sources(parallel PRIVATE
    "OpenMP.cpp"
)

set(PARALLEL_HEADER
    "MPI_environment.h"
    "OpenMP_buffer.h"
    "OpenMP.h")

target_include_directories(parallel PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

# target_compile_options(parallel PRIVATE -fPIC)

message("Required package of OpenMP is ${OpenMP_CXX_FOUND}")

# resource: https://gitlab.kitware.com/cmake/cmake/-/issues/21818
if(OpenMP_CXX_FOUND)
    target_link_libraries(parallel PUBLIC OpenMP::OpenMP_CXX)
endif()

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(parallel 
    INCLUDE_SUBDIR parallel 
    INCLUDES ${PARALLEL_HEADER})
    