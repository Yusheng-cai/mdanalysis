#include "SimulationState.h"

Calculation::Calculation(const CalculationInput& input)
:simstate_(input.simstate_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, vectorOutputs_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, vectorOutputNames_);

    ASSERT((vectorOutputs_.size() == vectorOutputNames_.size()), "The size of the vector outputs does not agree with the \
    size of names.");

    pack_.ReadVectorString("perIteroutputs", ParameterPack::KeyType::Optional, perIteroutputs_);
    pack_.ReadVectorString("perIteroutputNames", ParameterPack::KeyType::Optional, perIteroutputNames_);
    pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    for (int i=0;i<perIteroutputs_.size();i++)
    {
        ofsVector_.push_back(ofsptr(new std::ofstream));
    }

    for (int i=0;i<perIteroutputs_.size();i++)
    {
        ofsVector_[i]->open(perIteroutputNames_[i]);
    }
}

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

void Calculation::initializeResidueGroup(const std::string& residueName)
{
    ASSERT((! residueName.empty()), "The residue name is not provided.");

    // add the residue group
    addResidueGroup(residueName);

    // initialize the COM indice
    auto& res = getResidueGroup(residueName).getResidues();
    COMIndices_.resize(res[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1);

    // COMIndices are in 1-based counting 
    pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);
    for (int i=0;i<COMIndices_.size();i++)
    {
        COMIndices_[i] -= 1;
    }

    // resize the COM 
    COM_.resize(res.size());
}

void Calculation::printOutput()
{
    for (int i=0;i<vectorOutputs_.size();i++)
    {
        getOutputByName(vectorOutputs_[i])(vectorOutputNames_[i]);
    }

    closeAllOutputPerIter();
}

void Calculation::registerOutputFunction(std::string name, outputFunc func)
{
    auto it = MapNameToOutputFunction_.find(name);

    ASSERT((it == MapNameToOutputFunction_.end()), "The output with name " << name << " is already registered in calculation.");

    MapNameToOutputFunction_.insert(std::make_pair(name, func));
}

void Calculation::registerPerIterOutputFunction(std::string name, perIteroutputFunc func)
{
    auto it = MapNameToPerIterOutput_.find(name);

    ASSERT((it == MapNameToPerIterOutput_.end()), "The per iter output with name " << name << " is already registered in calculation.");

    MapNameToPerIterOutput_.insert(std::make_pair(name, func));
}

Calculation::perIteroutputFunc& Calculation::getIterOutputByName(std::string name)
{
    auto it = MapNameToPerIterOutput_.find(name);

    ASSERT((it != MapNameToPerIterOutput_.end()), "The per iter output with name " << name << " is not registered.");

    return it -> second;
}


Calculation::outputFunc& Calculation::getOutputByName(std::string name)
{
    auto it = MapNameToOutputFunction_.find(name);

    ASSERT((it != MapNameToOutputFunction_.end()), "The output with name " << name << " is not found within calculation");

    return it -> second;
}

void Calculation::printOutputOnStep()
{
    for (int i=0;i<perIteroutputs_.size();i++)
    {
        std::string name = perIteroutputs_[i];
        getIterOutputByName(name)(*ofsVector_[i]);
    }
}

void Calculation::closeAllOutputPerIter()
{
    for (int i=0;i<ofsVector_.size();i++)
    {
        ofsVector_[i]->close();
    }
}