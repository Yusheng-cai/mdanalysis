#include "liquid_crystal.h"

liquid_crystal::liquid_crystal(const OrderParametersInput& input)
:OrderParameters(input)
{
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);

    addAtomGroup(headgroupname_);
    addAtomGroup(tailgroupname_);

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    tailgroupsize_ = tailAG.getAtomGroupIndices().size();
    headgroupsize_ = headAG.getAtomGroupIndices().size();

    ASSERT((tailgroupsize_ == headgroupsize_), "The number of atoms passed in for headgroup is less than that of tail group."); 
}