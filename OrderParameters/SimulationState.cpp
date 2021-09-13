#include "SimulationState.h"

void SimulationState::registerAtomGroup(std::string name, AtomGroup& ag)
{
    //MapName2AtomGroup_.insert(std::make_pair(name,std::move(ag)));
    MapName2AtomGroup_.emplace(name, std::move(ag));
}

const AtomGroup& SimulationState::getAtomGroup(std::string name) const
{
    auto it = MapName2AtomGroup_.find(name);

    ASSERT((it != MapName2AtomGroup_.end()), "The name of AtomGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

AtomGroup& SimulationState::getAtomGroup(std::string name)
{
    auto it = MapName2AtomGroup_.find(name);

    ASSERT((it != MapName2AtomGroup_.end()), "The name of AtomGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

void SimulationState::registerResidueGroup(std::string name, ResidueGroup& res)
{
    MapName2ResidueGroup_.emplace(name, std::move(res));
}

const ResidueGroup& SimulationState::getResidueGroup(std::string name) const
{
    auto it = MapName2ResidueGroup_.find(name);

    ASSERT((it != MapName2ResidueGroup_.end()), "The name of ResidueGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

ResidueGroup& SimulationState::getResidueGroup(std::string name)
{
    auto it = MapName2ResidueGroup_.find(name);

    ASSERT((it != MapName2ResidueGroup_.end()), "The name of ResidueGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

void SimulationState::registerProbeVolume(const std::string name, const ProbeVolume* pv_ptr)
{
    auto it = MapName2PV.find(name);

    ASSERT((it == MapName2PV.end()), "The probe volume with name " << name << " is already registered.");

    MapName2PV.insert(std::make_pair(name, ProbeVolumePtr(const_cast<ProbeVolume*>(pv_ptr))));
}

const ProbeVolume& SimulationState::getProbeVolume(const std::string name) const
{
    auto it  = MapName2PV.find(name);

    ASSERT((it != MapName2PV.end()), "The probe volume with name " << name << " is not registered.");

    return *(it -> second);
}

ProbeVolume& SimulationState::accessProbeVolume(const std::string name)
{
    auto it = MapName2PV.find(name);

    ASSERT((it != MapName2PV.end()), "The probe volume with name " << name << " is not registered.");

    return *(it -> second);
}