#include "SimulationState.h"

void SimulationState::registerAtomGroup(std::string name, AtomGroup& ag)
{
    //MapName2AtomGroup_.insert(std::make_pair(name,std::move(ag)));
    MapName2AtomGroup_.emplace(name, std::move(ag));
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

ProbeVolume& SimulationState::getProbeVolume(const std::string name)
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

void SimulationState::registerCalculation(std::string name, const Calculation& calc)
{
    auto it = MapName2Calculation_.find(name);

    ASSERT((it == MapName2Calculation_.end()), "The calculation with name " << name << " is registered twice.");

    MapName2Calculation_.insert(std::make_pair(name, const_cast<Calculation*>(&calc)));
}

const Calculation& SimulationState::getCalculation(std::string name)
{

}