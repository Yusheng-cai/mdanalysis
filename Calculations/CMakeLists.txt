add_library(mdaCalculation SHARED)

target_sources(mdaCalculation PRIVATE
    "Bin.cpp"
    "Calculation.cpp"
    "Timeseries.cpp"
    "DataFileParser.cpp"
)

set(mdaCalculation_HEADER
    "Bin.h"
    "Calculation.h"
    "Timeseries.h"
    "DataFileParser.h"
)

target_include_directories(mdaCalculation
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(mdaCalculation PUBLIC mdatools mdaxdr mdaparallel)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaCalculation 
    INCLUDE_SUBDIR mdaCalculation
    INCLUDES ${mdaCalculation_HEADER})