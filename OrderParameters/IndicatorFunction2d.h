#pragma once
#include "tools/CommonTypes.h"
#include "IndicatorFunction.h"
#include "tools/Constants.h"

#include <vector>
#include <cmath>

class IndicatorFunction2d:public IndicatorFunction
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        IndicatorFunction2d(){};
        IndicatorFunction2d(Real sigma, Real ac, Real min, Real max);
        virtual ~IndicatorFunction2d(){};

        virtual void calculate(const Real& x, Real& h_x, Real& htilde_x, Real& dhtilde_dx) override;
        virtual void setLimits() override;

 
    private:
        Real min_;
        Real max_;

        Real lowlow_;
        Real lowhigh_;
        Real highlow_;
        Real highhigh_;
};