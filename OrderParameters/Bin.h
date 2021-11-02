#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>
#include <cmath>

class Bin
{
    public:
        using Real  = CommonTypes::Real;
        using Range = CommonTypes::Real2;

        Bin(const ParameterPack& pack);
        Bin(Real min, Real max, int numbins);

        int findBin(Real data);
        bool isInRange(Real data);

        Range getRange() const {return range_;}
        int getNumbins() const {return numbins_;}
        Real getStep() const {return step_;}
        Real getCenterLocationOfBin(int binNum) const;
        Real getLeftLocationOfBin(int binNum) const;
        Real getRightLocationOfBin(int binNum) const;


    private:
        Range range_;
        int numbins_;
        Real step_;
};