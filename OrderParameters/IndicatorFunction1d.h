#pragma once
#include "IndicatorFunction.h"

class IndicatorFunction1d:public IndicatorFunction
{
    public:
        IndicatorFunction1d(){};
        IndicatorFunction1d(Real sigma, Real ac, Real max);
        virtual ~IndicatorFunction1d(){};

        virtual void calculate(const Real& x, Real& h_x, Real& htilde_x, Real& dhtilde_dx) override;
        virtual void setLimits() override;
    
    private:
        Real max_;
        Real highhigh_;
        Real highlow_;
};