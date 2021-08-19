add_library(xdrC SHARED)

target_sources(xdrC PRIVATE 
                "trr_seek.c"
                "xdrfile_trr.c"
                "xdrfile_xtc.c"
                "xdrfile.c"
                "xtc_seek.c"
)

set(XDRC_HEADER
    "trr_seek.h"
    "xdrfile_trr.h"
    "xdrfile_xtc.h"
    "xdrfile.h"
    "xtc_seek.h"
)

# Add source files
target_include_directories(xdrC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(xdrC 
    INCLUDE_SUBDIR xdrC
    INCLUDES ${XDRC_HEADER})