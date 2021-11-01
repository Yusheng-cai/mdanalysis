#include "findNearStructure.h"

namespace CalculationRegistry
{
    registry_<findNearStructure> registerFnearStructure("findnearstructure");
}

findNearStructure::findNearStructure(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("mode", ParameterPack::KeyType::Optional, mode_);

    if (mode_ == "residue")
    {
        pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
        initializeResidueGroup(resName_);

        auto& res = getResidueGroup(resName_);
        COM_.resize(res.getResidues().size());
    }
    else if (mode_ == "atom")
    {
        pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, agName_);
        addAtomgroup(agName_);
    }
    else
    {
        ASSERT((true==false), "mode " << mode_ << " is not supported currently");
    }

    pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvName_);

    auto& pv = simstate_.getProbeVolume(pvName_);

    // register the per iter outputs
    registerPerIterOutputFunction("atomindices", [this](std::ofstream& ofs) -> void {this -> printAtomIndicesPerIter(ofs);});
}

void findNearStructure::calculateRes()
{
    AtomIndices_.clear();

    auto& pv = simstate_.getProbeVolume(pvName_);
    auto& res= simstate_.getResidueGroup(resName_).getResidues();

    for (int i=0;i<res.size();i++)
    {
        auto& r = res[i];

        Real3 C = CalculationTools::getCOM(r, simstate_, COMIndices_);
        COM_[i] = C;
    }

    for (int i=0;i<res.size();i++)
    {
        auto output = pv.calculate(COM_[i]);

        if (output.hx_ == 1)
        {
            for (int j=0;j<res[i].atoms_.size();j++)
            {
                AtomIndices_.push_back(res[i].atoms_[j].atomNumber_ -1);
            }
        }
    }
}

void findNearStructure::calculateAtom()
{
    AtomIndices_.clear();

    auto& pv  = simstate_.getProbeVolume(pvName_);
    auto& ag  = simstate_.getAtomGroup(agName_).getAtoms();

    for (int i=0;i<ag.size();i++)
    {
        auto output = pv.calculate(ag[i].position);

        if (output.hx_ == 1)
        {
            int Aindices = ag[i].index;
            AtomIndices_.push_back(Aindices);
        }
    }
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
    if (mode_ == "atom")
    {
        calculateAtom();
    }
    else if (mode_ == "residue")
    {
        calculateRes();
    }
}