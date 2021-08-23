add_library(mdaparallel SHARED)

target_sources(mdaparallel PRIVATE
    "OpenMP.cpp"
)

set(mdaPARALLEL_HEADER
    "MPI_environment.h"
    "OpenMP_buffer.h"
    "OpenMP.h")

target_include_directories(mdaparallel PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

# target_compile_options(mdaparallel PRIVATE -fPIC)

message("Required package of OpenMP is ${OpenMP_CXX_FOUND}")

# resource: https://gitlab.kitware.com/cmake/cmake/-/issues/21818
if(OpenMP_CXX_FOUND)
    target_link_libraries(mdaparallel PUBLIC OpenMP::OpenMP_CXX)
endif()

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaparallel 
    INCLUDE_SUBDIR mdaparallel 
    INCLUDES ${mdaPARALLEL_HEADER})
    