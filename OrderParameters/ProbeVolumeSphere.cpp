#include "ProbeVolumeSphere.h"

namespace ProbeVolumes
{
    static const registry_<ProbeVolumeSphere> registerSphere("sphere");
}

ProbeVolumeSphere::ProbeVolumeSphere(ProbeVolumeInput& input)
:ProbeVolume(input)
{
    input.ParamPack.ReadNumber("radius", ParameterPack::KeyType::Required, r_);
    input.ParamPack.ReadNumber("alpha_c", ParameterPack::KeyType::Optional, ac_);
    input.ParamPack.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);

    if (isDynamic())
    {
        addDynamicResidueGroup(dynamicResGroup_);
    }
    else
    {
        input.ParamPack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
    }

    func_ = IndicatorFunction1d(sigma_, ac_, r_); 

    setGeometry();
}

void ProbeVolumeSphere::update()
{
    if (isDynamic())
    {
        // obtain the residue group
        auto& res = getDynamicResidueGroup(dynamicResGroup_).getTotalResidue();

        std::vector<int> indices(res.atoms_.size());
        std::iota(indices.begin(), indices.end(), 0);

        // get the COM of the entire dynamic residue group
        Real3 COM = CalculationTools::getCOM(res, simstate_, indices);

        // update the center
        center_ = COM;

        #ifdef MY_DEBUG
        std::cout << "COM updated = " << center_[0] << " " << center_[1] << " " << center_[2] << std::endl;
        #endif 
    }
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