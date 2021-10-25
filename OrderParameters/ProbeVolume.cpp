#include "SimulationState.h"

ProbeVolume::ProbeVolume(ProbeVolumeInput& input):simstate_(input.simstate), simbox_(simstate_.getSimulationBox())
{
    // read in dynamic group
    isDynamic_ = input.ParamPack.ReadString("dynamicgroup", ParameterPack::KeyType::Optional, dynamicAtomName_);

    if (isDynamic())
    {
        addDynamicAtomGroup(dynamicAtomName_);
    }

    // read in alphaC and sigma
    input.ParamPack.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);
    input.ParamPack.ReadNumber("alpha_c", ParameterPack::KeyType::Optional, ac_);
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

void ProbeVolume::addDynamicAtomGroup(std::string name)
{
    auto it = MapAtomNameToIndex_.find(name);

    ASSERT((it == MapAtomNameToIndex_.end()), "The atom group with name " << name << " is registered twice.");

    auto& ag = simstate_.getAtomGroup(name);
    int size = AtomGroups_.size();
    MapAtomNameToIndex_.insert(std::make_pair(name, size));

    AtomGroups_.push_back(&ag);
}

AtomGroup& ProbeVolume::getDynamicAtomGroup(std::string name)
{
    auto it = MapAtomNameToIndex_.find(name);

    ASSERT((it != MapAtomNameToIndex_.end()), "The atomgroup with name " << name << " is not registered.");
    int index = it -> second;

    return *AtomGroups_[index];
}