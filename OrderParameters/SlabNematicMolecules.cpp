#include "SlabNematicMolecules.hpp"
#include "SimulationState.h"

namespace CalculationRegistry
{
    static registry_<SlabNematicMolecules> registerSlabNematicMolecules("SlabNematicMolecules");
}

SlabNematicMolecules::SlabNematicMolecules(const CalculationInput& input)
:Calculation(input)
{
   // initialize the bins first --> zbin and thetabin
    auto zbinPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);

    // default to be 30
    if(zbinPack != nullptr){
        zBin_    = Binptr(new Bin(*zbinPack));
        numzbins_= zBin_->getNumbins();
    }
    else
    {
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numzbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Required, above_);
        zBin_    = Binptr(new Bin());
        usingMinMax_ = true;
    }

    ReadInputs();

    pack_.ReadNumber("left_max", ParameterPack::KeyType::Optional, left_max_);
    pack_.ReadNumber("right_max", ParameterPack::KeyType::Optional, right_max_);
    pack_.ReadNumber("right_alphaC", ParameterPack::KeyType::Optional, right_alphaC_);
    pack_.ReadNumber("left_alphaC", ParameterPack::KeyType::Optional, left_alphaC_);

    func_left_ = IndicatorFunction1d(left_sigma_, left_alphaC_, left_max_);
    func_right_ = IndicatorFunction1d(right_sigma_, right_alphaC_, right_max_);

    vector_stilde_.resize(numzbins_, 0.0);
    ResidueLocationPerBin_.resize(numzbins_,0.0);
    numResiduePerBin_.resize(numzbins_, 0.0);

    registerOutputFunction("slab_stilde", [this](std::string name) -> void {this -> printAverageVectorStilde(name);});
}


void SlabNematicMolecules::ReadInputs()
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional,headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional,tailIndex_);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    initializeResidueGroup(residueGroupName_);

    // correct for the indexing on head and tail index 
    headIndex_--;
    tailIndex_--;

    // find the direction
    directionIndex_ = MapdirectionToIndex_.find(direction_) -> second;
}


void SlabNematicMolecules::binUsingMinMax(){
    Real slight_shift=1e-3;
    std::vector<Real> ZdirectionNum;

    for (int i=0;i<COM_.size();i++){
        Real val = COM_[i][directionIndex_];

        if (val > above_){
            ZdirectionNum.push_back(val);
        }
    }

    auto maxit = std::max_element(ZdirectionNum.begin(), ZdirectionNum.end());
    auto minit = std::min_element(ZdirectionNum.begin(), ZdirectionNum.end());

    Real max = *maxit + slight_shift;
    Real min = *minit - slight_shift;
    Range2 range = {{min, max}};

    std::cout << "Max = " << max << " Min = " << min << std::endl;

    zBin_ -> update(range, numzbins_);

    for (int i=0;i<numzbins_;i++){
        ResidueLocationPerBin_[i] += zBin_->getCenterLocationOfBin(i);
    }

}

void SlabNematicMolecules::calculate(){
   // obtain the residue group
    auto& res = getResidueGroup(residueGroupName_).getResidues();
    vector_stilde_iter_.clear();
    vector_stilde_iter_.resize(numzbins_,0.0);
    uij_.clear();
    uij_.resize(res.size());
    COM_.clear();
    COM_.resize(res.size());


    // find all the COM of the residues in the system
    #pragma omp parallel for 
    for (int i=0;i<res.size();i++){
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

        Real3 distance;
        Real dist_sq;
        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, dist_sq);
        LinAlg3x3::normalize(distance);
        uij_[i] = distance;
    }

    // bin using min max of the molecules if needed 
    if (usingMinMax_){
        binUsingMinMax();
    }

    for (int i=0;i<res.size();i++){
        Real dotProduct = LinAlg3x3::DotProduct(uij_[i], arr_);
        Real s_left_i, stilde_left_i, dstilde_left_i_dx;
        Real s_right_i, stilde_right_i, dstilde_right_i_dx;

        func_left_.calculate(dotProduct, \
                                                s_left_i, stilde_left_i, dstilde_left_i_dx);
        func_right_.calculate(dotProduct, \
                                s_right_i, stilde_right_i, dstilde_right_i_dx);
        if (zBin_->isInRange(COM_[i][directionIndex_])){
            int idx = zBin_->findBin(COM_[i][directionIndex_]);
            numResiduePerBin_[idx] += 1;
            if ((1-s_right_i) == 1 || s_left_i == 1){
                vector_stilde_[idx] += 1;
                vector_stilde_iter_[idx] += 1;
            }
        }
    }
}

void SlabNematicMolecules::finishCalculate(){
    int num_frames = simstate_.getTotalFrames();

    vector_stilde_ = vector_stilde_ / num_frames;
    ResidueLocationPerBin_ = ResidueLocationPerBin_ / num_frames;
    numResiduePerBin_ = numResiduePerBin_ / num_frames;
}

void SlabNematicMolecules::printAverageVectorStilde(std::string inputfname){
    std::ofstream ofs;
    ofs.open(inputfname);

    for (int i=0;i<zBin_->getNumbins();i++){
        ofs << i << " " << ResidueLocationPerBin_[i] << " " << numResiduePerBin_[i] << " " \
        << vector_stilde_[i] << "\n";
    }

    ofs.close();
}


