#include "SimulationState.h"

ProbeVolume::ProbeVolume(ProbeVolumeInput& input):simstate_(input.simstate), simbox_(simstate_.getSimulationBox())
{
    // read in dynamic group
    isDynamic_ = input.ParamPack.ReadString("dynamicgroup", ParameterPack::KeyType::Optional, dynamicResGroup_);
};

void ProbeVolume::addDynamicResidueGroup(std::string name)
{
    auto it = MapResidueNameToIndex_.find(name);

    ASSERT((it == MapResidueNameToIndex_.end()), "The residue with name " << name << " is registered twice.");

    auto& res = simstate_.getResidueGroup(name);
    int size = ResidueGroups_.size();
    MapResidueNameToIndex_.insert(std::make_pair(name, size));

    ResidueGroups_.push_back(&res);
}

ResidueGroup& ProbeVolume::getDynamicResidueGroup(std::string name)
{
    auto it = MapResidueNameToIndex_.find(name);

    ASSERT((it != MapResidueNameToIndex_.end()), "The residue with name " << name << " is not registered.");
    int index = it -> second;

    return *ResidueGroups_[index];
}