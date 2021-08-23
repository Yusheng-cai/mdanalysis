#include "ProbeVolumeCylinder.h"

ProbeVolumeCylinder::ProbeVolumeCylinder(ProbeVolumeInput& input)
:ProbeVolume(input)
{

}

void ProbeVolumeCylinder::setGeometry(Real rmax, Real zmax, Real ac, Real sigma)
{
    rrange_[0] = 0;
    rrange_[1] = rmax;

    zrange_[0] = -zmax/2;
    zrange_[1] = zmax/2;

    ac_ = ac;
    sigma_ = sigma;

    zfunc_ = IndicatorFunction2d(sigma_, ac_, zrange_[0], zrange_[1]);
    rfunc_ = IndicatorFunction1d(sigma_, ac_, rmax);
}

ProbeVolumeOutput ProbeVolumeCylinder::calculate(const Real3& x) const
{
    Real3 diff;

    for (int i=0;i<3;i++)
    {
        diff[i] = x[i] - center_[i];
    }

    // radius is x^2 + y^2
    Real radius = std::sqrt(diff[1]*diff[1] + diff[0]*diff[0]);

    Real dhtildedr_;
    Real htilder_;
    Real hr_;

    // calculate the indicator function in the r direction
    rfunc_.calculate(radius,hr_, htilder_,dhtildedr_);


}