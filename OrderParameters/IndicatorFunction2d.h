#include "tools/CommonTypes.h"
#include "IndicatorFunction.h"

#include <vector>
#include <cmath>

class IndicatorFunction2d:public IndicatorFunction
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        IndicatorFunction2d();
        IndicatorFunction2d(Real sigma, Real ac);
    
    private:
        Real sigma_=0.01;
        Real ac_=0.02;

        Real k_;
        Real k1_;
        Real k2_;
};