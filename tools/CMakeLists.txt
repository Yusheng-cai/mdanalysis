add_library(tools SHARED)

target_sources(tools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "Frame.cpp"
    "XdrWrapper.cpp"
    "TrrFile.cpp"
)

set(TOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
)

add_subdirectory(libxdr)

target_include_directories(tools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(tools PRIVATE xdr)
# target_compile_options(tools PRIVATE -fPIC)


#include(InstallLibraryWithStandardSetup)
#InstallLibraryWithStandardSetup(tools 
#    INCLUDE_SUBDIR tools
#    INCLUDES ${TOOLS_HEADER})

