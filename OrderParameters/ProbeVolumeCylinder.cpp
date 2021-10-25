#include "ProbeVolumeCylinder.h"

ProbeVolumeCylinder::ProbeVolumeCylinder(ProbeVolumeInput& input)
:ProbeVolume(input)
{
    //ASSERT((! isDynamic()), "The cylinder probe volume does not yet support dynamic atom group.");
}

void ProbeVolumeCylinder::setGeometry(Real rmax, Real zmax, Real ac, Real sigma)
{
    rrange_[0] = 0;
    rrange_[1] = rmax;

    // z goes from 0 to zmax
    zrange_[0] = 0;
    zrange_[1] = zmax;

    ac_ = ac;
    sigma_ = sigma;

    zfunc_ = IndicatorFunction2d(sigma_, ac_, zrange_[0], zrange_[1]);
    rfunc_ = IndicatorFunction1d(sigma_, ac_, rmax);
}

ProbeVolumeOutput ProbeVolumeCylinder::calculate(const Real3& x) const
{
    // radius is x^2 + y^2
    Real radius = std::sqrt(x[1]*x[1] + x[2]*x[2]);

    Real dhtildedr_;
    Real htilder_;
    Real hr_;

    // calculate the indicator function in the r direction
    rfunc_.calculate(radius,hr_, htilder_,dhtildedr_);

    Real dhtildedz_;
    Real htildez_;
    Real hz_;
    zfunc_.calculate(x[0],hz_, htildez_, dhtildedz_);

    // Calculate the non-coarse grained indicator for cylinder
    Real hcylinder = hz_*hr_;
    Real htildecylinder = htildez_*htilder_;

    // Calculate derivative of the htilde with respect to each of the coordinate
    Real3 dhtildedc_;

    // if radius is larger than almost 0, we calculate derivatives normally
    // if radius is 0, then the derivative is 0
    if (radius >= 1e-7)
    {
        // dhtildedx
        dhtildedc_[0] = dhtildedr_/radius*x[0]*htildez_;

        // dhtildedy 
        dhtildedc_[1] = dhtildedr_/radius*x[1]*htildez_;

        // dhtildedz
        dhtildedc_[2] = dhtildedz_*htilder_;
    }
    else
    {
        dhtildedc_[0] = 0.0;
        dhtildedc_[1] = 0.0;

        dhtildedc_[2] = dhtildedz_*htilder_;
    }

    ProbeVolumeOutput output;

    output.hx_ = hcylinder;
    output.htilde_x_ = htildecylinder;
    output.dhtilde_dx_ = dhtildedc_;

    return output;
}