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

    Real3 getCOM(const Molecule::residue& residues, const SimulationState& simstate, std::vector<int>& indices_);
};