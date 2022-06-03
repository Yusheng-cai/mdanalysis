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

        /*
        Constructor for the Bin object 
        There are 3 types of constructors

        1. It can take in a ParameterPack object, which means basically a block in the input file
        2. It can directly take in the range [min, max] and the number of bins between the min and max
        3. Empty constructor

        */
        Bin(const ParameterPack& pack);
        Bin(Range range, int numbins);
        Bin(){};

        int findBin(Real data);
        bool isInRange(Real data);
        void update(Range range, int numbins);

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