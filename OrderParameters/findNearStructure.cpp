#include "findNearStructure.h"

namespace CalculationRegistry
{
    registry_<findNearStructure> registerFnearStructure("findnearstructure");
}

findNearStructure::findNearStructure(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvName_);

    auto& pv = simstate_.getProbeVolume(pvName_);
    ASSERT((pv.isDynamic()), "The probe volume must be dynamic for this calculation.");

    // register the per iter outputs
    registerPerIterOutputFunction("atomindices", [this](std::ofstream& ofs) -> void {this -> printAtomIndicesPerIter(ofs);});

    // add the residue groups
    initializeResidueGroup(resname_);

    auto& res = getResidueGroup(resname_).getResidues();
    COM_.resize(res.size());
}

void findNearStructure::printAtomIndicesPerIter(std::ofstream& ofs)
{
    ofs << simstate_.getStep() << " ";
    for (int i=0;i<AtomIndices_.size();i++)
    {
        ofs << AtomIndices_[i] << " ";
    }
    ofs << "\n";
}

void findNearStructure::calculate()
{
    AtomIndices_.clear();

    auto& res = getResidueGroup(resname_).getResidues();
    auto& pv  = simstate_.getProbeVolume(pvName_);

    for (int i=0;i<res.size();i++)
    {
        Real3 com = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = com;
    }

    for (int i=0;i<COM_.size();i++)
    {
        auto output = pv.calculate(COM_[i]);

        if (output.hx_ == 1)
        {
            for (int j=0;j<res[i].atoms_.size();j++)
            {
                int Aindices = res[i].atoms_[j].atomNumber_ - 1;
                AtomIndices_.push_back(Aindices);
            }
        }
    }
}