add_library(tools SHARED)

add_subdirectory(libxdr)

target_sources(tools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "XdrWrapper.cpp"
)

set(TOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "XdrWrapper.h"
)

target_include_directories(tools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(tools PUBLIC xdr)

# include(InstallLibraryWithStandardSetup)
# InstallLibraryWithStandardSetup(tools 
#     INCLUDE_SUBDIR tools
#     INCLUDES ${TOOLS_HEADER})

