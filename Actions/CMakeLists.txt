add_library(mdaActions SHARED)

target_sources(mdaActions PRIVATE
    "MorsePotential.cpp"
    "mda_actions.cpp"
    "NanoparticleGeneration.cpp"
    "LennardJones.cpp"
)

set(mdaActions_HEADER
    "MorsePotential.hpp"
    "mda_actions.hpp"
    "NanoparticleGeneration.hpp"
    "LennardJones.hpp"
)


target_include_directories(mdaActions
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(mdaActions PUBLIC mdaOrderParameters mdatools mdaxdr mdaparallel ${FFTW3_LIBRARIES} ${Boost_LIBRARIES})

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdaActions 
    INCLUDE_SUBDIR mdaActions
    INCLUDES ${mdaActions_HEADER})
