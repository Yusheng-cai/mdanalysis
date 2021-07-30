#include "IndicatorFunction2d.h"

IndicatorFunction2d::IndicatorFunction2d(Real sigma, Real ac, Real min, Real max)
:min_(min), max_(max),IndicatorFunction(sigma, ac)
{
    setLimits();
}


void IndicatorFunction2d::calculate(const Real& x, Real& h_x, Real& htilde_x, Real& dhtilde_dx) const
{
    if (x >= min_ &&  x<= max_)
    {
        h_x = 1;
    }
    else
    {
        h_x = 0;
    }

    // within the inner region
    if (x > lowhigh_ && x<highlow_)
    {
        htilde_x = 1;
        dhtilde_dx = 0;
    }

    if (x < lowlow_ || x > highhigh_)
    {
        htilde_x = 0;
        dhtilde_dx = 0;
    }

    if (x >= lowlow_ && x <= lowhigh_)
    {
        // The 1 from the last theta function is already added
        htilde_x = k1_*std::erf((x - min_)/(std::sqrt(2)*sigma_)) - k2_*(x - min_) + 0.5;

        Real tilde1 = truncatedGaussian(max_ - x);
        Real tilde2 = truncatedGaussian(min_ - x);

        dhtilde_dx = -(tilde1 - tilde2);
    }

    if (x >= highlow_ && x <= highhigh_)
    {
        // The 1 from the last theta funciton is already added
        htilde_x = k1_*std::erf((max_ - x)/(std::sqrt(2)*sigma_)) - k2_*(max_ - x) + 0.5;

        Real tilde1 = truncatedGaussian(max_ - x);
        Real tilde2 = truncatedGaussian(min_ - x);

        dhtilde_dx = -(tilde1 - tilde2);
    }
}

void IndicatorFunction2d::setLimits()
{
    lowlow_ = min_ - ac_;
    lowhigh_= min_ + ac_;
    highlow_= max_ - ac_;
    highhigh_= max_ + ac_;
}