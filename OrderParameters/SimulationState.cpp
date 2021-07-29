#include "SimulationState.h"

void SimulationState::registerAtomGroup(std::string name, AtomGroup& ag)
{
    MapName2AtomGroup_.insert(std::make_pair(name,ag));
}

const AtomGroup& SimulationState::getAtomGroup(std::string name) const
{
    auto it = MapName2AtomGroup_.find(name);

    ASSERT((it != MapName2AtomGroup_.end()), "The name of AtomGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}