#include "IndicatorFunction.h"

IndicatorFunction::IndicatorFunction(Real sigma, Real ac)
:sigma_(sigma), ac_(ac)
{
    sigma2_ = sigma_*sigma_;
    ac2_ = ac_*ac_;

    calculateFactors();
}

IndicatorFunction::Real IndicatorFunction::truncatedGaussian(Real alpha) const
{
    if ((ac_ - std::abs(alpha))>= 0)
    {
        Real ret = 1.0/k_*(std::exp(-alpha*alpha/(2*sigma2_)) - std::exp(-ac2_/(2*sigma2_)));
        return ret;
    }
    else
    {
        return 0.0;
    }
}

void IndicatorFunction::calculateFactors()
{
    k_ = std::sqrt(2*Constants::PI*sigma2_)*std::erf(ac_/std::sqrt(2*sigma2_)) - 2*ac_*std::exp(-ac2_/(2*sigma2_));
    k1_ = 1.0/k_*std::sqrt(Constants::PI*sigma2_/2);
    k2_ = 1.0/k_*std::exp(-ac2_/(2*sigma2_));
}