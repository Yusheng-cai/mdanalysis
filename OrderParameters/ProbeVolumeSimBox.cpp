#include "ProbeVolumeSimBox.h"

namespace ProbeVolumes
{
    static const registry_<ProbeVolumeSimBox> register_Simbox("simbox");
}

ProbeVolumeSimBox::ProbeVolumeSimBox(ProbeVolumeInput& input)
:ProbeVolume(input)
{
    ASSERT((! isDynamic()), "The simulation box probe volume does not support dynamic atom group.");
}

ProbeVolumeOutput ProbeVolumeSimBox::calculate(const Real3& x) const
{
    ProbeVolumeOutput output;

    output.dhtilde_dx_ = {{0,0,0}};
    output.htilde_x_ = 1;
    output.hx_ = 1;

    return output;
}