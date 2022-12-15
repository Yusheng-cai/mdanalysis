#pragma once

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
#include <cstdio>

namespace mda_tools{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    void WriteGroFile(const std::string& filename, const std::vector<Real3>& vecpos, const std::vector<std::string>& aname, const std::vector<std::string>& resname, \
                      const Real3& boxsize);
    void readGroFile(const std::string& filename, std::vector<Real3>& vecpos, std::vector<std::string>& aname, std::vector<std::string>& resname);

    void readFileLines(const std::string& filename, std::vector<std::string>& lines);
}