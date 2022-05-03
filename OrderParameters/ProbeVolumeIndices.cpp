#include "ProbeVolumeIndices.h"

namespace CalculationRegistry
{
    registry_<ProbeVolumeIndices> registerPVIndices("ProbeVolumeIndices");
}

ProbeVolumeIndices::ProbeVolumeIndices(const CalculationInput& input)
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
    else if (mode_ == "atomres")
    {
        pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
        addResidueGroup(resName_);
    }
    else
    {
        ASSERT((true==false), "mode " << mode_ << " is not supported currently");
    }

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // register the per iter outputs
    registerPerIterOutputFunction("atomindices", [this](std::ofstream& ofs) -> void {this -> printAtomIndicesPerIter(ofs);});
}

void ProbeVolumeIndices::calculateRes()
{
    AtomIndices_.clear();

    auto& res= simstate_.getResidueGroup(resName_).getResidues();

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
        bool inpv = false;
        for (auto pv: NotInprobevolumes_)
        {
            auto output = pv -> calculate(COM_[i]);

            if(output.hx_ == 1)
            {
                inpv = true;
            }
        }

        if (! inpv)
        {
            // TODO: this does not work if we have union
            for (auto pv : probevolumes_)
            {
                auto output = pv->calculate(COM_[i]);
                if (output.hx_ == 1)
                {
                    insideIndices.push_back(i);
                }
            }
        }
    }

    for (int i=0;i<insideIndices.size();i++)
    {
        int index = insideIndices[i];
        for (int j=0;j<res[index].atoms_.size();j++)
        {
            AtomIndices_.push_back(res[index].atoms_[j].atomNumber_-1);
        }
    }
}

void ProbeVolumeIndices::calculateAtom()
{
    AtomIndices_.clear();

    auto& ag  = simstate_.getAtomGroup(agName_).getAtoms();

    // find the inside indices 
    std::vector<int> insideIndices;
    for (int i=0;i<ag.size();i++)
    {
        bool inpv = false;
        for (auto pv: NotInprobevolumes_)
        {
            auto output = pv -> calculate(ag[i].position);

            if(output.hx_ == 1)
            {
                inpv = true;
            }
        }

        if (! inpv)
        {
            // TODO: this does not work if we have union
            for (auto pv : probevolumes_)
            {
                auto output = pv->calculate(ag[i].position);
                if (output.hx_ == 1)
                {
                    AtomIndices_.push_back(ag[i].index);
                }
            }
        }
    }
}

void ProbeVolumeIndices::calculateAtomRes()
{
    AtomIndices_.clear();

    auto& res = simstate_.getResidueGroup(resName_).getResidues();

    // find the inside indices 
    std::vector<int> ResidueIndices;
    std::map<int,bool> MapresIndTobool;
    for (int i=0;i<res.size();i++)
    {
        bool inside =false;
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            if (inside)
            {
                break;
            }

            auto& atom = res[i].atoms_[j];
            bool inpv = false;
            for (auto pv: NotInprobevolumes_)
            {
                auto output = pv -> calculate(atom.positions_);

                if(output.hx_ == 1)
                {
                    inpv = true;
                }
            }

            if (! inpv)
            {
                // TODO: this does not work if we have union
                for (auto pv : probevolumes_)
                {
                    auto output = pv->calculate(atom.positions_);
                    if (output.hx_ == 1)
                    {
                        ResidueIndices.push_back(i);
                        inside = true;
                        break;
                    }
                }
            }
        }
    }

    // let's check for duplicates in residue indices 
    std::map<int, bool> DuplicateMap;
    for (int i=0;i<ResidueIndices.size();i++)
    {
        auto it = DuplicateMap.find(ResidueIndices[i]);

        ASSERT((it == DuplicateMap.end()), "There is duplicate residue indices, something wrong with the code.");
    }

    // Then let's add the atom indices 
    for (int i=0;i<ResidueIndices.size();i++)
    {
        int index = ResidueIndices[i];
        for (int j=0;j<res[index].atoms_.size();j++)
        {
            AtomIndices_.push_back(res[index].atoms_[j].atomNumber_-1);
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
    else if (mode_ == "atomres")
    {
        calculateAtomRes();
    }
}