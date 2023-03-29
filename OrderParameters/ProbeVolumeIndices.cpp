#include "ProbeVolumeIndices.h"

namespace CalculationRegistry
{
    registry_<ProbeVolumeIndices> registerPVIndices("ProbeVolumeIndices");
}

ProbeVolumeIndices::ProbeVolumeIndices(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("mode", ParameterPack::KeyType::Optional, mode_);

    if (mode_ == "residue"){
        pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
        initializeResidueGroup(resName_);

        auto& res = getResidueGroup(resName_);
        COM_.resize(res.getResidues().size());
    }
    else if (mode_ == "atom"){
        pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, agName_);
        addAtomgroup(agName_);
    }
    else{
        ASSERT((true==false), "mode " << mode_ << " is not supported currently");
    }

    // initialize probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // register the per iter outputs
    registerPerIterOutputFunction("atomindices", [this](std::ofstream& ofs) -> void {this -> printAtomIndicesPerIter(ofs);});
}

void ProbeVolumeIndices::calculateRes()
{
    // clear atom indices 
    AtomIndices_.clear();
    auto& res= simstate_.getResidueGroup(resName_).getResidues();

    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        auto& r = res[i];

        Real3 C = CalculationTools::getCOM(r, simstate_, COMIndices_);
        COM_[i] = C;
    }

    // find the inside indices 
    std::vector<int> insideIndices;
    for (int i=0;i<COM_.size();i++)
    {
        if (isInPV(COM_[i]))
        {
            for (int j=0;j<res[i].atoms_.size();j++)
            {
                AtomIndices_.push_back(res[i].atoms_[j].atomNumber_);
            }
        }
    }
}

void ProbeVolumeIndices::calculateAtom()
{
    AtomIndices_.clear();

    auto& ag  = simstate_.getAtomGroup(agName_).accessAtoms();

    // find the inside indices 
    std::vector<int> insideIndices;
    for (int i=0;i<ag.size();i++)
    {
        if (isInPV(ag[i].position))
        {
            AtomIndices_.push_back(ag[i].index);
        }
    }
}

void ProbeVolumeIndices::printAtomIndicesPerIter(std::ofstream& ofs)
{
    ofs << simstate_.getStep() << " ";
    for (int i=0;i<AtomIndices_.size();i++)
    {
        ofs << AtomIndices_[i] << " ";
    }
    ofs << "\n";
}

void ProbeVolumeIndices::calculate()
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