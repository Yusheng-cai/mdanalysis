#include "SlabResidueProperty.h"

namespace CalculationRegistry
{
    registry_<SlabResidueProperty> registerSlabResidueProperty("SlabResidueProperty");
}

SlabResidueProperty::SlabResidueProperty(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    pack_.ReadNumber("direction", ParameterPack::KeyType::Optional, direction_);

    initializeResidueGroup(resname_);

    auto binPack = pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);
    if(binPack != nullptr){
        bin_    = Binptr(new Bin(*binPack));
        numbins_= bin_->getNumbins();
        dz_     = bin_->getStep();
    }
    else{
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Required, above_);
        bin_    = Binptr(new Bin());
        usingMinMax_ = true;
        dz_     = 0.0;
    }
    ResidueLocationPerBin_.resize(numbins_,0.0);
    histogram_.resize(numbins_,0.0);

    registerOutputFunction("histogram", [this](std::string name) -> void {this -> printHistogram(name);});
}

void SlabResidueProperty::binUsingMinMax(const std::vector<Real3>& positions)
{
    Real slight_shift=1e-3;
    std::vector<Real> ZdirectionNum;

    for (int i=0;i<positions.size();i++){
        Real val = positions[i][direction_];

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

    bin_ -> update(range, numbins_);
    dz_ = bin_->getStep();

    for (int i=0;i<numbins_;i++){
        ResidueLocationPerBin_[i] += bin_->getCenterLocationOfBin(i);
    }
}

void SlabResidueProperty::printHistogram(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# binLocation   num\n";
    
    for (int i=0;i<numbins_;i++){
        ofs << ResidueLocationPerBin_[i] << " " << histogram_[i] << "\n";
    }

    ofs.close();
}

void SlabResidueProperty::calculate(){
    // first calculate the COM 
    auto res = getResidueGroup(resname_);
    for (int i=0;i<COM_.size();i++){
        COM_[i] = calcCOM(res[i]);
    }

    if (usingMinMax_){
        binUsingMinMax(COM_);
    }

    for (int i=0;i<COM_.size();i++){
        if (bin_->isInRange(COM_[i][direction_])){
            int index = bin_ -> findBin(COM_[i][direction_]);
            histogram_[index] += 1;
        }
    }
}

void SlabResidueProperty::finishCalculate(){
    int numframes = simstate_.getTotalFrames();

    histogram_ = histogram_ / numframes;

    if (usingMinMax_){
        ResidueLocationPerBin_ = ResidueLocationPerBin_ / numframes;
    }
    else{
        for (int i=0;i<numbins_;i++){
            ResidueLocationPerBin_[i] = bin_->getCenterLocationOfBin(i);
        }
    }
}