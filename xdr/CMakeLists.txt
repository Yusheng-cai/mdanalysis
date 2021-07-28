add_subdirectory(libxdr)
add_library(xdr SHARED)

target_sources(xdr PRIVATE
    "XdrWrapper.cpp"
    "TrrFile.cpp"
    "Frame.cpp"
    "XtcFile.cpp"
)

target_include_directories(xdr
                        PUBLIC ${CMAKE_SOURCE_DIR})


set(XDR_HEADER
    "Frame.h"
    "TrrFile.h"
    "XdrWrapper.h"
    "XtcFile.h"
)

target_link_libraries(xdr PRIVATE tools xdrC)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(xdr 
    INCLUDE_SUBDIR xdr
    INCLUDES ${XDR_HEADER})