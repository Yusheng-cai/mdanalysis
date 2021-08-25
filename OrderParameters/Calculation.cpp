#include "Calculation.h"

Calculation::Calculation(const CalculationInput& input)
:simstate_(input.simstate_)
{}

void Calculation::addAtomgroup(std::string name)
{
    // Check if the AtomGroupName is not in the map
    auto it = MapAtomGroupNameToIndex_.find(name);
    ASSERT((it == MapAtomGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the AtomGroup reference from simulation state
    AtomGroup& ag = simstate_.getAtomGroup(name);

    int index = AtomGroups_.size();
    AtomGroups_.push_back(&ag);

    MapAtomGroupNameToIndex_.insert(std::make_pair(name, index));
}

void Calculation::addResidueGroup(std::string name)
{
    // Check if the ResidueGroupName is not in the map
    auto it = MapResidueGroupNameToIndex_.find(name);
    ASSERT((it == MapResidueGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the ResidueGroup reference from simulation state
    ResidueGroup& rg = simstate_.getResidueGroup(name);

    int index = ResidueGroups_.size();
    ResidueGroups_.push_back(&rg);

    MapResidueGroupNameToIndex_.insert(std::make_pair(name, index));
}

const ResidueGroup& Calculation::getResidueGroup(std::string name) const 
{
    auto it  = MapResidueGroupNameToIndex_.find(name);

    ASSERT((it != MapResidueGroupNameToIndex_.end()), "The residue with name " << name << " is not registered.");

    return *ResidueGroups_[it->second]; 
}