#pragma once
#include "xdr/GroFile.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "xdr/MoleculeStructs.h"
#include "SimulationState.h"
#include "parallel/OpenMP_buffer.h"

#include <vector>
#include <cmath>
#include <array>

namespace CalculationTools
{
    using Real  = CommonTypes::Real;
    using Real3 = CommonTypes::Real3;

    std::vector<Real3> getCOM(const std::vector<Molecule::residue>& residues, const SimulationState& simstate);
};