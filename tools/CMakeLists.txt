add_library(mdatools SHARED)

target_sources(mdatools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "FileSystem.cpp"
    "Algorithm.cpp"
)

set(mdaTOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "Registry.h"
    "FileSystem.h"
    "Algorithm.h"
)

target_include_directories(mdatools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)
# target_compile_options(mdatools PRIVATE -fPIC)


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(mdatools 
    INCLUDE_SUBDIR mdatools
    INCLUDES ${mdaTOOLS_HEADER})

