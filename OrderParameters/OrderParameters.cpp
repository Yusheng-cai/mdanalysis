#include "OrderParameters.h"

OrderParameters::OrderParameters(const OrderParametersInput& input)
:simstate_(input.simstate), simbox_(input.simstate.getSimulationBox())
{
    input.pack_.ReadString("name",ParameterPack::KeyType::Required,name_);
}

void OrderParameters::registerOutput(std::string name, OutputValue::ValueFunction func)
{
    OutputValue output(name, func);

    OutputValues.insert(std::make_pair(name, output)); 
}

const OutputValue& OrderParameters::getOutput(std::string name) const
{
    auto it = OutputValues.find(name);

    ASSERT((it != OutputValues.end()), "The name " << name << " is not available.");

    return it -> second;
}

int OrderParameters::getAtomGroupIndex(std::string& name) const
{
    auto it = MapAtomGroupNameToIndex_.find(name);

    ASSERT(( it != MapAtomGroupNameToIndex_.end() ), "The Atom Group with name " << name << " does not exist.");

    int index = it->second;

    return index;
}

void OrderParameters::addAtomGroup(std::string name)
{
    // Check if the AtomGroupName is not in the map
    auto it = MapAtomGroupNameToIndex_.find(name);
    ASSERT((it == MapAtomGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the AtomGroup reference from simulation state
    AtomGroup& ag = simstate_.getAtomGroup(name);

    int index = OPAtomGroups_.size();
    OPAtomGroups_.push_back(&ag);

    MapAtomGroupNameToIndex_.insert(std::make_pair(name, index));

    derivativeOutputs_.push_back(AtomDerivativesSet());

    ASSERT((derivativeOutputs_.size() == OPAtomGroups_.size()), "The derivative size does not match with number of atomgroups.");
}

const AtomDerivativesSet& OrderParameters::getDerivatives(std::string& name) const
{
    int index = getAtomGroupIndex(name); 

    return derivativeOutputs_[index];
}

AtomDerivativesSet& OrderParameters::accessDerivatives(std::string& name)
{
    int index = getAtomGroupIndex(name); 

    return derivativeOutputs_[index];
}

const AtomGroup& OrderParameters::getAtomGroup(std::string& name) const
{
    int index = getAtomGroupIndex(name);    

    return *OPAtomGroups_[index];
}

AtomGroup& OrderParameters::accessAtomGroup(std::string& name)
{
    int index = getAtomGroupIndex(name);

    return const_cast<AtomGroup&>(*OPAtomGroups_[index]);
}

void OrderParameters::clearDerivative(std::string& name)
{
    accessDerivatives(name).clear();
}

void OrderParameters::clearDerivativesOutputs()
{
    for (int i=0;i<derivativeOutputs_.size();i++)
    {
        derivativeOutputs_[i].clear();

        derivativeOutputs_[i].setMasterObject();
    }
}