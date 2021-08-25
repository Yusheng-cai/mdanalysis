#include "QtensorZ.h"

QtensorZ::QtensorZ(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residuegroup", ParameterPack::KeyType::Required, residueName_);
    input.pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    input.pack_.ReadArrayNumber("range", ParameterPack::KeyType::Required, range_);

    auto it = MapNameToDirection.find(direction_);
    ASSERT((it != MapNameToDirection.end()), "The direction " << direction_ << " is not recognized.");
    index_ = it -> second;

    // add the residue group to the system
    addResidueGroup(residueName_);
}

void QtensorZ::calculate()
{
    // obtain the residue group by its name
    const auto& res = getResidueGroup(residueName_);

    // obtain the COM
    std::vector<Real3> COM = CalculationTools::getCOM(res.getResidues(),simstate_);
}