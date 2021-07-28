add_library(xdrC SHARED)

target_sources(xdrC PRIVATE 
                "trr_seek.c"
                "xdrfile_trr.c"
                "xdrfile_xtc.c"
                "xdrfile.c"
                "xtc_seek.c"
)

# Add source files
target_include_directories(xdrC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})