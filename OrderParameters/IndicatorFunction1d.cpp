#include "IndicatorFunction1d.h"

IndicatorFunction1d::IndicatorFunction1d(Real sigma, Real ac, Real max)
:max_(max),IndicatorFunction(sigma, ac)
{
    setLimits();
}

void IndicatorFunction1d::calculate(const Real& x, Real& h_x, Real& htilde_x, Real& dhtilde_dx) const
{
    if (x <= max_)
    {
        h_x = 1;
    }
    else
    {
        h_x = 0;
    }

    if (x > highhigh_)
    {
        htilde_x = 0;
        dhtilde_dx = 0;
    }

    if (x < highlow_)
    {
        htilde_x = 1;
        dhtilde_dx = 0;
    }

    if (x >= highlow_ && x<= highhigh_)
    {
        htilde_x = k1_*std::erf((max_ - x)/(std::sqrt(2)*sigma_)) - k2_*(max_ - x) + 0.5; 

        dhtilde_dx = - truncatedGaussian(max_ - x);
    }
}

void IndicatorFunction1d::setLimits()
{
    highlow_ = max_ - ac_;
    highhigh_ = max_ + ac_;
}