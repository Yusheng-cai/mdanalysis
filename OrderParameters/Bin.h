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

        int findBin(Real data);
        bool isInRange(Real data);

        Range getRange() const {return range_;}
        int getNumbins() const {return numbins_;}
        Real getStep() const {return step_;}
        int getDimension() const {return dimension_;}
        Real getLocationOfBin(int binNum) const;


    private:
        std::vector<Real> bins_;

        Range range_;
        int numbins_;
        int dimension_ = 1;
        Real step_;
};