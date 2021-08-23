#include "Bin.h"

Bin::Bin(const ParameterPack& pack)
{
    pack.ReadArrayNumber("range", ParameterPack::KeyType::Required, range_);
    pack.ReadNumber("numbins", ParameterPack::KeyType::Required,numbins_);
    pack.ReadNumber("dimension", ParameterPack::KeyType::Optional, dimension_);

    step_ = (range_[1] - range_[0])/numbins_;
}

int Bin::findBin(Real x)
{
    int index = std::floor((x - range_[0])/step_);

    return index;
}

bool Bin::isInRange(Real data)
{
    if (data >= range_[0] && data <= range_[1])
    {
        return true;
    }
    else
    {
        return false;
    }
}

Bin::Real Bin::getLocationOfBin(int binNum) const 
{
    return range_[0] + (binNum + 0.5) * step_;
}
