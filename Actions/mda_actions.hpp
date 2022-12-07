#include "tools/CommandLineArguments.h"
#include "tools/CommonOperations.h"
#include "tools/CommonOperations.h"
#include "tools/InputParser.h"
#include "NanoparticleGeneration.hpp"
#include "xdr/GroFile.h"

#include <string>
#include <iostream>
#include <vector>
#include <array>

namespace mda_actions{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    void generateNP(CommandLineArguments& cmd);

    void FindSurfaceSiO2(CommandLineArguments& cmd);
}