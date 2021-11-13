#include "cost.h"

namespace CalculationRegistry
{
    static registry_<Cost> registerPcost("Cost");
}

Cost::Cost(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,tailIndex_);

    input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    input.pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    input.pack_.ReadString("probevolume", ParameterPack::KeyType::Required, ProbeVolumeName_);
    input.pack_.ReadNumber("numbins", ParameterPack::KeyType::Optional, numBins_);

    // initialize the bin pointer 
    Range range_ = {{0.0, 1.0}};
    Bin_ = Binptr(new Bin(range_,numBins_));

    // resize the histogram 
    histogram_.resize(numBins_, 0.0);

    auto bfPack = input.pack_.findParamPack("BetaFactors", ParameterPack::KeyType::Optional);

    if (bfPack != nullptr)
    {
        bf_ = BFptr(new BetaFactorWriter(const_cast<ParameterPack&>(*bfPack)));
    }

    headIndex_--;
    tailIndex_--;

    // add the residue group
    addResidueGroup(residueGroupName_);

    auto& res = getResidueGroup(residueGroupName_).getResidues();
    COMIndices_.resize(res[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1);
    input.pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);

    for (int i=0;i<COMIndices_.size();i++)
    {
        COMIndices_[i] -= 1;
    }

    COM_.resize(res.size());

    // resize the molecular director size
    uij_.resize(res.size());

    registerOutputFunction("histogram", [this](std::string name) -> void {this -> printhistogram(name);});
    registerPerIterOutputFunction("costheta", [this](std::ofstream& ofs) -> void {this -> printavgCosthetaPerIter(ofs);});
    registerPerIterOutputFunction("Ntilde", [this](std::ofstream& ofs) -> void { this -> printNtilde(ofs);});

    // register for output file
    registerOutputFileOutputs("costheta", [this](void) -> Real {return this -> getCostheta2();});
}


void Cost::calculate()
{
    BetaFactors_.clear();
    
    auto& res = getResidueGroup(residueGroupName_).getResidues();
    auto& pv  = simstate_.getProbeVolume(ProbeVolumeName_);

    int totalatomNumbers = getResidueGroup(residueGroupName_).getAtomSize();
    BetaFactors_.resize(totalatomNumbers);
    std::fill(BetaFactors_.begin(), BetaFactors_.end(), -1.0);


    // find all the COM of the residues in the system
    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        Real3 com = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = com;

        Real3 distance_;
        Real dist_sq;
        Real3 headPos_ = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos_ = res[i].atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(headPos_, tailPos_, distance_, dist_sq);

        Real3 normalized_dir = Qtensor::normalize_director(distance_);
        uij_[i] = normalized_dir;
    }

    // check which COM are inside the probevolume
    std::vector<int> InsideIndices;
    for (int i=0;i<COM_.size();i++)
    {
        auto pvOutput = pv.calculate(COM_[i]);
        if (pvOutput.hx_ == 1)
        {
            InsideIndices.push_back(i);
        }
    }

    avgCostheta_ = 0.0;

    for (int i=0;i<InsideIndices.size();i++)
    {
        int k = InsideIndices[i];
        Real cost = Qtensor::vec_dot(uij_[k], arr_);
        Real cost2 = cost * cost;

        // average the cos squared theta 
        avgCostheta_ += cost2;

        ASSERT((cost >= -1 && cost <= 1), "cosine(theta) is not within range of -1 and 1");

        for (int j=0;j<res[k].atoms_.size();j++)
        {
            int atomNum = res[k].atoms_[j].atomNumber_;

            BetaFactors_[atomNum-1] = std::pow(cost,2.0);
        }

        // start binning 
        ASSERT((Bin_->isInRange(cost2)), "cosine(theta)^2 is not within range of 0 and 1.");
        int binNum = Bin_->findBin(cost2);
        histogram_[binNum] += 1;
    }

    avgCostheta_ /= InsideIndices.size();
}

void Cost::printavgCosthetaPerIter(std::ofstream& ofs)
{
    // print time and cos theta 
    Real framenum = simstate_.getFrameNumber();

    ofs << framenum << " " << avgCostheta_ << "\n";
}

void Cost::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= numFrames;
    }
}

void Cost::printNtilde(std::ofstream& ofs)
{
    // obtain the residue 
    auto& res = getResidueGroup(residueGroupName_).getTotalResidue();
    auto& pv  = simstate_.getProbeVolume(ProbeVolumeName_);

    Real Ntilde = 0.0;

    #pragma omp parallel
    {
        Real Ntilde_local = 0.0;
        #pragma omp for
        for (int i=0;i<res.atoms_.size();i++)
        {
            auto pvOutput = pv.calculate(res.atoms_[i].positions_);

            Ntilde_local += pvOutput.htilde_x_;
        }

        #pragma omp critical
        {
            Ntilde += Ntilde_local;
        }
    }

    int numframe = simstate_.getFrameNumber();

    ofs << numframe << " " << Ntilde << " " << avgCostheta_ << std::endl;
}

void Cost::printOutputOnStep()
{
    if (bf_.get() != nullptr)
    {
        bf_ -> write(simstate_.getFrameNumber(), BetaFactors_);
    }
}

void Cost::printhistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# index location value\n";

    for (int i=0;i<histogram_.size();i++)
    {
        Real location = Bin_ -> getCenterLocationOfBin(i);
        ofs_ << i << " " << location << " " << histogram_[i] << "\n";
    }

    ofs_.close();
}