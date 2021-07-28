#include "ProbeVolumeBox.h"

namespace ProbeVolumeRegistry
{
    static const registry_<ProbeVolumeBox> register_box("box");
}

ProbeVolumeBox::ProbeVolumeBox(ProbeVolumeInput& input)
:ProbeVolume(input), state(input.simstate)
{
    input.ParamPack.ReadArrayNumber("xrange", ParameterPack::KeyType::Required, xrange_);
    input.ParamPack.ReadArrayNumber("yrange", ParameterPack::KeyType::Required, yrange_);
    input.ParamPack.ReadArrayNumber("zrange", ParameterPack::KeyType::Required, zrange_);
    input.ParamPack.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);
    input.ParamPack.ReadNumber("alphac", ParameterPack::KeyType::Optional, ac_);

    // Initialize the indicator functions
    dx = xrange_[1] - xrange_[0];
    dy = yrange_[1] - yrange_[0];
    dz = zrange_[1] - zrange_[0];

    center_[0] = 0.5*(xrange_[0] + xrange_[1]);
    center_[1] = 0.5*(yrange_[0] + yrange_[1]);
    center_[2] = 0.5*(zrange_[0] + zrange_[1]);

    func_[0] = IndicatorFunction2d(sigma_, ac_, -dx/2, dx/2);
    func_[1] = IndicatorFunction2d(sigma_, ac_, -dy/2, dy/2);
    func_[2] = IndicatorFunction2d(sigma_, ac_, -dz/2, dz/2);
}

ProbeVolumeOutput ProbeVolumeBox::calculate(const Real3& x)
{
}
