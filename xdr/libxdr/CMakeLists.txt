add_library(mdaxdrC SHARED)

target_sources(mdaxdrC PRIVATE 
                "trr_seek.c"
                "xdrfile_trr.c"
                "xdrfile_xtc.c"
                "xdrfile.c"
                "xtc_seek.c"
)

set(mdaXDRC_HEADER
    "trr_seek.h"
    "xdrfile_trr.h"
    "xdrfile_xtc.h"
    "xdrfile.h"
    "xtc_seek.h"
)

# Add source files
target_include_directories(mdaxdrC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaxdrC 
    INCLUDE_SUBDIR mdaxdrC
    INCLUDES ${mdaXDRC_HEADER})