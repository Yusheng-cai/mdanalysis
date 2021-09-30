#include "Bin.h"

Bin::Bin(const ParameterPack& pack)
{
    pack.ReadArrayNumber("range", ParameterPack::KeyType::Required, range_);
    pack.ReadNumber("numbins", ParameterPack::KeyType::Required,numbins_);

    step_ = (range_[1] - range_[0])/numbins_;
}

int Bin::findBin(Real x)
{
    int index;

    Real epsilon = 1e-5;
    if ( std::abs(x - range_[0]) <= epsilon)
    {
        index = 0;
    }
    else
    {
        index = std::ceil((x-range_[0])/step_) - 1;
    }

    ASSERT((index >= 0 && index < numbins_), "Bin index out of range, it is equal to " << index);

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

Bin::Real Bin::getCenterLocationOfBin(int binNum) const 
{
    return range_[0] + (binNum + 0.5) * step_;
}

Bin::Real Bin::getLeftLocationOfBin(int binNum) const 
{
    return range_[0] + binNum * step_;
}

Bin::Real Bin::getRightLocationOfBin(int binNum) const 
{
    return range_[0] + (binNum + 1) * step_; 
}
