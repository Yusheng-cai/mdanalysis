#include "ProbeVolumeSphere.h"

namespace ProbeVolumes
{
    static const registry_<ProbeVolumeSphere> registerSphere("sphere");
}

ProbeVolumeSphere::ProbeVolumeSphere(ProbeVolumeInput& input)
:ProbeVolume(input)
{
    input.ParamPack.ReadNumber("radius", ParameterPack::KeyType::Required, r_);
    input.ParamPack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
    input.ParamPack.ReadNumber("alpha_c", ParameterPack::KeyType::Optional, ac_);
    input.ParamPack.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);

    func_ = IndicatorFunction1d(sigma_, ac_, r_); 

    setGeometry();
}

void ProbeVolumeSphere::setGeometry()
{
    rmax_ = r_ + ac_;
}

ProbeVolumeOutput ProbeVolumeSphere::calculate(const Real3& x) const
{
    ProbeVolumeOutput output;

    Real3 dist;
    Real dist_sq;
    Real r;

    simbox_.calculateDistance(x, center_, dist, dist_sq);

    r = std::sqrt(dist_sq);

    Real h_r;
    Real htilde_r;
    Real dhtilde_dr;
    Real3 dhtilde;
    func_.calculate(r, h_r, htilde_r, dhtilde_dr);

    // The derivative is dh/dx = dh/dr*dr/dx
    // r = sqrt(x^2 + y^2 + z^2), dr/dx = x/r
    // dh/dx = x/r * dh/dr

    for (int i=0;i<3;i++)
    {
        dhtilde[i] = dhtilde_dr/r*dist[i];
    }

    output.htilde_x_ = htilde_r;
    output.dhtilde_dx_ = dhtilde;
    output.hx_ = h_r;

    return output;
}