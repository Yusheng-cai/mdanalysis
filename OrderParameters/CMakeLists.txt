add_library(mdaOrderParameters SHARED)

target_sources(mdaOrderParameters PRIVATE
    "IndicatorFunction.cpp"
    "IndicatorFunction1d.cpp"
    "IndicatorFunction2d.cpp"
    "OrderParameters.cpp"
    "ProbeVolume.cpp"
    "ProbeVolumeBox.cpp"
    "ProbeVolumeSphere.cpp"
    "ProbeVolumeTiltedCylinder.cpp"
    "ProbeVolumeCylinder.cpp"
    "ProbeVolumeSimBox.cpp"
    "SimulationBox.cpp"
    "SimulationState.cpp"
    "AtomGroup.cpp"
    "Output_values.cpp"
    "Indus.cpp"
    "Driver.cpp"
    "Output_files.cpp"
    "LiquidCrystal.cpp"
    "P2.cpp"
    "P2tilde.cpp"
    "P2cos.cpp"
    "AtomGroupParsingStrategy.cpp"
    "AtomDerivativesSet.cpp"
    "ResidueGroup.cpp"
    "Bin.cpp"
    "Calculation.cpp"
    "DataFileParser.cpp"
    "SlabQtensor.cpp"
    "Timeseries.cpp"
    "RDFresidue.cpp"
    "SlabOrientation.cpp"
    "cost.cpp"
    "BetaFactorWriter.cpp"
    "ProbeVolumeSimBox.cpp"
    "Dipole.cpp"
    "LinAlgTools.cpp"
    "ProbeVolumeIndices.cpp"
    "FFT.cpp"
    "MSD.cpp"
    "LocalOrder.cpp"
    "gcost.cpp"
    "grcost.cpp"
    "grcostorientation.cpp"
    "DistributionGivenDimer.cpp"
    "SRE.cpp"
    "CellGrid.cpp"
    "SRE_dimer.cpp"
    "OrientationalDistribution.cpp"
    "PositionalOrientationalAvg.cpp"
)

set(mdaOrderParameters_HEADER
    "IndicatorFunction.h"
    "IndicatorFunction1d.h"
    "IndicatorFunction2d.h"
    "OrderParameters.h"
    "ProbeVolume.h"
    "ProbeVolumeBox.h"
    "ProbeVolumeSphere.h"
    "ProbeVolumeCylinder.h"
    "ProbeVolumeTiltedCylinder.h"
    "ProbeVolumeSimBox.h"
    "SimulationState.h"
    "AtomGroup.h"
    "Output_values.h"
    "Indus.h"
    "Driver.h"
    "Output_files.h"
    "LiquidCrystal.h"
    "P2.h"
    "P2tilde.h"
    "P2cos.h"
    "AtomGroupParsingStrategy.h"
    "AtomDerivativesSet.h"
    "ResidueGroup.h"
    "Bin.h"
    "Calculation.h"
    "DataFileParser.h"
    "SlabQtensor.h"
    "Timeseries.h"
    "RDFresidue.h"
    "SlabOrientation.h"
    "cost.h"
    "BetaFactorWriter.h"
    "ProbeVolumeSimBox.h"
    "Dipole.h"
    "LinAlgTools.h"
    "FFT.h"
    "MSD.h"
    "LocalOrder.h"
    "gcost.h"
    "grcost.h"
    "grcostorientation.h"
    "DistributionGivenDimer.h"
    "SRE.h"
    "CellGrid.h"
    "SRE_dimer.h"
    "OrientationalDistribution.h"
    "ProbeVolumeIndices.h"
    "PositionalOrientationalAvg.h"
)


target_include_directories(mdaOrderParameters
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_include_directories(mdaOrderParameters PUBLIC ${FFTW3_INCLUDE_DIRS})

message("In src, ${FFTW3_LIBRARIES}")

target_link_libraries(mdaOrderParameters PUBLIC mdatools mdaxdr mdaparallel ${FFTW3_LIBRARIES})

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaOrderParameters 
    INCLUDE_SUBDIR mdaOrderParameters
    INCLUDES ${mdaOrderParameters_HEADER})
