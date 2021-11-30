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
    input.pack_.ReadNumber("numbins", ParameterPack::KeyType::Optional, numBins_);

    // normalize the vector 
    LinAlg3x3::normalize(arr_);

    // initialize the bin pointer 
    Range range_ = {{0.0, 1.0}};
    Bin_ = Binptr(new Bin(range_,numBins_));

    // resize the histogram 
    histogram_.resize(numBins_, 0.0);
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

    // initialize the probe volumes 
    initializeProbeVolumes();
    initializeNotInProbeVolumes();
}


void Cost::calculate()
{
    auto& res = getResidueGroup(residueGroupName_).getResidues();
    std::vector<ProbeVolume*> probevolumes;
    for (int i=0;i<probevolumeNames_.size();i++)
    {
        auto& pv  = simstate_.getProbeVolume(probevolumeNames_[i]);
        probevolumes.push_back(&pv);
    }

    // resize the histogram
    histogram_.resize(numBins_, 0.0);

    int totalatomNumbers = getResidueGroup(residueGroupName_).getAtomSize();

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
        for (auto pv : probevolumes)
        {
            auto pvOutput = pv->calculate(COM_[i]);
            if (pvOutput.hx_ == 1)
            {
                InsideIndices.push_back(i);
                break;
            }
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

        if ( !(cost >= -1 && cost <= 1))
        {
            std::cout << "WARNING cosine(theta) is not within range of -1 and 1 and it is " << cost << std::endl;
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