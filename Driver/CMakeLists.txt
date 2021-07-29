add_library(Driver SHARED)

target_sources(Driver PRIVATE
    "Driver.cpp"
)

set(DRIVER_HEADER
    "Driver.h"
)

target_include_directories(Driver
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(Driver PUBLIC tools xdr OrderParameters)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(Driver
    INCLUDE_SUBDIR Driver
    INCLUDES ${DRIVER_HEADER})

