#include "NematicMolecules.hpp"
#include "SimulationState.h"

namespace CalculationRegistry
{
    registry_<NematicMolecules> registerNematicMolecules("NematicMolecules");
}

NematicMolecules::NematicMolecules(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, ResidueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,HeadIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,TailIndex_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    pack_.Readbool("usedirector", ParameterPack::KeyType::Optional, useDirector_);
    pack_.ReadNumber("left_max", ParameterPack::KeyType::Optional, left_max_);
    pack_.ReadNumber("right_max", ParameterPack::KeyType::Optional, right_max_);
    pack_.ReadNumber("right_alphaC", ParameterPack::KeyType::Optional, right_alphaC_);
    pack_.ReadNumber("left_alphaC", ParameterPack::KeyType::Optional, left_alphaC_);

    func_left_ = IndicatorFunction1d(left_sigma_, left_alphaC_, left_max_);
    func_right_ = IndicatorFunction1d(right_sigma_, right_alphaC_, right_max_);

    HeadIndex_--;
    TailIndex_--;

    // initialize Residue
    initializeResidueGroup(ResidueGroupName_);

    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    registerOutputFileOutputs("stilde", [this](void) -> Real {return this->get_stilde();});

}

void NematicMolecules::calculate(){
    auto& res    = getResidueGroup(ResidueGroupName_).getResidues();
    int AtomSize = simstate_.getTotalNumberAtoms();
    stilde_= 0.0;
    COM_.clear();
    COM_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());
    AtomIndices_.clear();

    std::vector<int> ResidueIndex_;

    #pragma omp parallel for 
    for (int i=0;i<res.size();i++){
        COM_[i] = calcCOM(res[i]);

        Real3 HeadPos = res[i].atoms_[HeadIndex_].positions_;
        Real3 TailPos = res[i].atoms_[TailIndex_].positions_;
        Real3 distance;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(HeadPos, TailPos,distance, dist_sq);
        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    for (int i=0;i<res.size();i++){
        Real dotProduct = LinAlg3x3::DotProduct(uij_[i], arr_);
        Real s_left_i, stilde_left_i, dstilde_left_i_dx;
        Real s_right_i, stilde_right_i, dstilde_right_i_dx;


        func_left_.calculate(dotProduct, \
                                                s_left_i, stilde_left_i, dstilde_left_i_dx);
        func_right_.calculate(dotProduct, \
                                s_right_i, stilde_right_i, dstilde_right_i_dx);

        stilde_ += (stilde_left_i + 1 - stilde_right_i);

        if (s_right_i == 1 || s_left_i == 1){
            ResidueIndex_.push_back(i);
        }
    }

    for (int i=0;i<ResidueIndex_.size();i++){
        int ind = ResidueIndex_[i];

        for (int j=0;j<res[i].atoms_.size();j++){
            AtomIndices_.push_back(res[i].atoms_[j].atomNumber_);
        }
    }
}