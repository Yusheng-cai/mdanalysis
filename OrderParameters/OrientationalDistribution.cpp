#include "OrientationalDistribution.h"

namespace CalculationRegistry
{
    registry_<OrientationalDistribution> registerOrientationalDistribution("OrientationalDistribution");
}

OrientationalDistribution::OrientationalDistribution(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, ResidueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,HeadIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,TailIndex_);
    pack_.ReadNumber("numbins", ParameterPack::KeyType::Required, NumBins_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    HeadIndex_--;
    TailIndex_--;

    // initialize Residue
    initializeResidueGroup(ResidueGroupName_);

    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // initialize the bins 
    CosThetaBin_ = Binptr(new Bin(CosThetaRange_, NumBins_));
    CosThetaSquaredBin_ = Binptr(new Bin(CosThetaSquaredRange_, NumBins_));
    PCosTheta_.resize(NumBins_,0.0);
    PCosThetaSquared_.resize(NumBins_,0.0);

    // register the output functions 
    registerOutputs();
    registerOutputfile();
}

void OrientationalDistribution::registerOutputfile()
{
    registerOutputFileOutputs("costheta", [this](void)-> Real {return this -> getAvgCostheta();});
    registerOutputFileOutputs("costhetasquared", [this](void)-> Real {return this -> getAvgCosthetasquared();});
}


void OrientationalDistribution::registerOutputs()
{
    registerOutputFunction("Distribution", [this](std::string name) -> void {this -> PrintDistribution(name);});
}

void OrientationalDistribution::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();
    Real CosThetaTot = 0.0;
    Real CosSquaredThetaTot = 0.0;

    for (int i=0;i<NumBins_;i++)
    {
        PCosTheta_[i] = PCosTheta_[i]/numFrames;
        PCosThetaSquared_[i] = PCosThetaSquared_[i] / numFrames;

        CosThetaTot += PCosTheta_[i];
        CosSquaredThetaTot += PCosThetaSquared_[i];
    }

    for (int i=0;i<NumBins_;i++)
    {
        PCosTheta_[i] = PCosTheta_[i] / CosThetaTot;
        PCosThetaSquared_[i] = PCosThetaSquared_[i] / CosSquaredThetaTot;
    }
}

void OrientationalDistribution::PrintDistribution(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Bin\tcostheta\tcos2theta\tP(costheta)\tP(cos2theta)\n";

    for (int i=0;i<NumBins_;i++)
    {
        ofs <<CosThetaBin_->getCenterLocationOfBin(i) << " " << CosThetaSquaredBin_->getCenterLocationOfBin(i) << " " \
        << PCosTheta_[i] << " " << PCosThetaSquared_[i] << "\n";
    }
    ofs.close();
}

void OrientationalDistribution::calculate()
{
    auto& res = getResidueGroup(ResidueGroupName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());
    AvgCostheta_ = 0.0;
    AvgCosthetasquared_= 0.0;

    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 HeadPos = res[i].atoms_[HeadIndex_].positions_;
        Real3 TailPos = res[i].atoms_[TailIndex_].positions_;
        Real3 distance;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(HeadPos, TailPos,distance, dist_sq);
        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    // find atoms within probe volume
    std::vector<int> AtomIndices = InsidePVIndices(COM_);

    // iterate over the atom indices
    for (int i=0;i<AtomIndices.size();i++)
    {
        int k = AtomIndices[i];
        Real cost = LinAlg3x3::DotProduct(uij_[k], arr_);
        Real costsq = cost * cost;

        AvgCostheta_ += cost;
        AvgCosthetasquared_ += costsq;

        int CosThetaBinNum = CosThetaBin_->findBin(cost);
        int CosThetaSquaredBinNum = CosThetaSquaredBin_->findBin(costsq);

        PCosTheta_[CosThetaBinNum] += 1;
        PCosThetaSquared_[CosThetaSquaredBinNum] += 1;
    }

    AvgCosthetasquared_ = AvgCosthetasquared_ * 1.0/AtomIndices.size();
    AvgCostheta_ = AvgCostheta_ * 1.0/AtomIndices.size();
}