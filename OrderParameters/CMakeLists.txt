add_library(OrderParameters SHARED)

target_sources(OrderParameters PRIVATE
    "IndicatorFunction.cpp"
    "IndicatorFunction1d.cpp"
    "IndicatorFunction2d.cpp"
    "OrderParameters.cpp"
    "ProbeVolume.cpp"
    "ProbeVolumeBox.cpp"
    "SimulationBox.cpp"
    "SimulationState.cpp"
    "AtomGroup.cpp"
    "Output_values.cpp"
    "Indus.cpp"
    "ProbeVolumeRegistry.cpp"
    "Driver.cpp"
    "Output_files.cpp"
)

set(OrderParameters_HEADER
    "IndicatorFunction.h"
    "IndicatorFunction1d.h"
    "IndicatorFunction2d.h"
    "OrderParameters.h"
    "ProbeVolume.h"
    "ProbeVolumeBox.h"
    "SimulationState.h"
    "AtomGroup.h"
    "Output_values.h"
    "Indus.h"
    "ProbeVolumeRegistry.h"
    "Driver.h"
    "Output_files.h"
)

target_include_directories(OrderParameters
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(OrderParameters PUBLIC tools xdr)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(OrderParameters 
    INCLUDE_SUBDIR OrderParameters
    INCLUDES ${OrderParameters_HEADER})