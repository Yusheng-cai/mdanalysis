add_subdirectory(libxdr)
add_library(mdaxdr SHARED)

target_sources(mdaxdr PRIVATE
    "XdrWrapper.cpp"
    "TrrFile.cpp"
    "Frame.cpp"
    "XtcFile.cpp"
    "GroFile.cpp"
    "TopologyReader.cpp"
)

target_include_directories(mdaxdr
                        PUBLIC ${CMAKE_SOURCE_DIR})


set(mdaXDR_HEADER
    "Frame.h"
    "TrrFile.h"
    "XdrWrapper.h"
    "XtcFile.h"
    "GroFile.h"
    "TopologyReader.h"
)

target_link_libraries(mdaxdr PRIVATE mdatools mdaxdrC)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaxdr 
    INCLUDE_SUBDIR mdaxdr
    INCLUDES ${mdaXDR_HEADER})