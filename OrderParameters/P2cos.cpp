#include "P2cos.h"

namespace OrderParametersRegistry
{
    registry_<P2cos> registerP2cos("p2cos");
}

P2cos::P2cos(const OrderParametersInput& input)
:OrderParameters(input)
{
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);
    input.pack_.ReadArrayNumber<Real,3>("director", ParameterPack::KeyType::Required, n_);

    addAtomGroup(headgroupname_);
    addAtomGroup(tailgroupname_);

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    tailgroupsize_ = tailAG.getAtomGroupIndices().size();
    headgroupsize_ = headAG.getAtomGroupIndices().size();

    ASSERT((tailgroupsize_ == headgroupsize_), "The number of atoms passed in for headgroup is less than that of tail group.");

    registerOutput("p2",[this](void)-> Real {return this->getP2cos();});
}

void P2cos::calculate()
{

}