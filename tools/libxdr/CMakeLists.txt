# Prefer to keep flags private)
add_library(xdr)

target_sources(xdr PRIVATE 
                "trr_seek.c"
                "xdrfile_trr.c"
                "xdrfile_xtc.c"
                "xdrfile.c"
                "xtc_seek.c"
)

# Add source files
target_include_directories(xdr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})