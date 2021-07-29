#include "ProbeVolumeRegistry.h"

void ProbeVolumeRegistry::registerProbeVolume(const std::string name, const ProbeVolume* pv_ptr)
{
    auto it = MapName2PV.find(name);

    ASSERT((it == MapName2PV.end()), "The probe volume with name " << name << " is already registered.");

    MapName2PV.insert(std::make_pair(name, ProbeVolumePtr(const_cast<ProbeVolume*>(pv_ptr))));
}

ProbeVolume& ProbeVolumeRegistry::getProbeVolume(const std::string name) const
{
    auto it  = MapName2PV.find(name);

    ASSERT((it != MapName2PV.end()), "The probe volume with name " << name << " is not registered.");

    return *(it -> second);
}