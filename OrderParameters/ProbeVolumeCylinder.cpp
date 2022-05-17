#include "ProbeVolumeCylinder.h"

ProbeVolumeCylinder::ProbeVolumeCylinder(ProbeVolumeInput& input)
:ProbeVolume(input)
{
}

void ProbeVolumeCylinder::setGeometry(Real rmin, Real rmax, Real zmax, Real ac, Real sigma)
{
    rrange_[0] = rmin;
    rrange_[1] = rmax;

    // z goes from 0 to zmax
    zrange_[0] = 0;
    zrange_[1] = zmax;

    ac_ = ac;
    sigma_ = sigma;

    zfunc_ = IndicatorFunction2d(sigma_, ac_, zrange_[0], zrange_[1]);
    rfunc_ = IndicatorFunction2d(sigma_, ac_, rrange_[0], rrange_[1]);
}

ProbeVolumeOutput ProbeVolumeCylinder::calculate(const Real3& x) const
{
    // radius is x^2 + y^2
    Real radius = std::sqrt(x[1]*x[1] + x[2]*x[2]);

    Real dhtildedr;
    Real htilder;
    Real hr;

    // calculate the indicator function in the r direction
    rfunc_.calculate(radius,hr, htilder,dhtildedr);

    Real dhtildedz;
    Real htildez;
    Real hz;
    zfunc_.calculate(x[0],hz, htildez, dhtildedz);

    // Calculate the non-coarse grained indicator for cylinder
    Real hcylinder = hz*hr;
    Real htildecylinder = htildez*htilder;

    // Calculate derivative of the htilde with respect to each of the coordinate
    Real3 dhtildedc;

    // if radius is larger than almost 0, we calculate derivatives normally
    // if radius is 0, then the derivative is 0
    if (radius >= 1e-7)
    {
        // dhtildedx
        dhtildedc[0] = dhtildedr/radius*x[0]*htildez;

        // dhtildedy 
        dhtildedc[1] = dhtildedr/radius*x[1]*htildez;

        // dhtildedz
        dhtildedc[2] = dhtildedz*htilder;
    }
    else
    {
        dhtildedc[0] = 0.0;
        dhtildedc[1] = 0.0;

        dhtildedc[2] = dhtildedz*htilder;
    }

    ProbeVolumeOutput output;

    output.hx_ = hcylinder;
    output.htilde_x_ = htildecylinder;
    output.dhtilde_dx_ = dhtildedc;

    return output;
}