#include "DistributionGivenDimer.h"


namespace CalculationRegistry
{
    registry_<DistributionGivenDimer> registerDistributionGivenDimer("DistGivenDimer");
}

DistributionGivenDimer::DistributionGivenDimer(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
    initializeResidueGroup(resName_);
    auto& res = getResidueGroup(resName_).getResidues();
    numres_ = res.size();
    AngleWithSurface_.resize(numres_,0.0);
    COM_.resize(numres_);
    uij_.resize(numres_);

    // read head index and tail index
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    // read in the binning information
    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_ -> getNumbins();
    histogram_.resize(numbins_,0.0);
    histogramNotDimer_.resize(numbins_,0.0);

    // read in the surface normal
    pack_.ReadArrayNumber("surfacenormal", ParameterPack::KeyType::Optional, surfaceNormal_);

    // read in the cosine theta max and r max 
    pack_.ReadNumber("cosmax", ParameterPack::KeyType::Required, cosmax_);
    pack_.ReadNumber("rmax", ParameterPack::KeyType::Required, rmax_);

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});
    registerOutputFunction("histogramnotdimer", [this](std::string name) -> void {this->printHistogramNotdimer(name);});
    registerPerIterOutputFunction("dimers", [this](std::ofstream& ofs) -> void {this -> printNumDimerPerIter(ofs);});
    registerOutputFileOutputs("dimers", [this](void) -> Real {return this->getNumDimers();});
    registerPerIterOutputFunction("dimerperres", [this](std::ofstream& ofs) -> void {this -> printNumDimerPerResiduePerIter(ofs);});
}

void DistributionGivenDimer::calculate()
{
    auto& res = getResidueGroup(resName_).getResidues();
    std::vector<int> InsideIndices;
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

        Real3 headpos = res[i].atoms_[headindex_].positions_;
        Real3 tailpos = res[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distsq);
        
        LinAlg3x3::normalize(distance);

        AngleWithSurface_[i] = LinAlg3x3::DotProduct(distance, surfaceNormal_);
        uij_[i] = distance;
    }

    // find the inside indices 
    InsideIndices = InsidePVIndices(COM_);
    int size = InsideIndices.size();
    std::vector<std::vector<Real>> pairDistances(size, std::vector<Real>(size,0.0));
    std::vector<std::vector<Real>> costhetaPair(size, std::vector<Real>(size,0.0));
    DimerPerResidue_.clear();
    DimerPerResidue_.resize(size,0.0);

    #pragma omp parallel for
    for (int i=0;i<size;i++)
    {
        for (int j=i+1;j<size;j++)
        {
            int index1 = InsideIndices[i];
            int index2 = InsideIndices[j];

            Real3 distance;
            Real distsq;

            simstate_.getSimulationBox().calculateDistance(COM_[index1], COM_[index2], distance, distsq);

            pairDistances[i][j] = std::sqrt(distsq);

            Real costheta = LinAlg3x3::DotProduct(uij_[index1], uij_[index2]);
            costhetaPair[i][j] = costheta;
        }
    }

    numDimersPerIter_ = 0;
    for (int i=0;i<size;i++)
    {
        for (int j=i+1;j<size;j++)
        {
            if (pairDistances[i][j] <= rmax_ && costhetaPair[i][j] <= cosmax_)
            {
                DimerPerResidue_[i] += 1;
                DimerPerResidue_[j] += 1;
                numDimersPerIter_ += 1;
            }
        }
    }

    for (int i=0;i<DimerPerResidue_.size();i++)
    {
        int index = InsideIndices[i];

        Real angle = AngleWithSurface_[index];

        if (bin_ -> isInRange(angle))
        {
            int binnum = bin_ -> findBin(angle);
            histogram_[binnum] += 1;

            if (DimerPerResidue_[i] > 0)
            {
                histogram_[index] += 1;
            }
            else
            {
                histogramNotDimer_[index] += 1;
            }
        }
    }
}

void DistributionGivenDimer::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= numFrames;
    }
}

void DistributionGivenDimer::printHistogramNotdimer(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<histogramNotDimer_.size();i++)
    {
        ofs << bin_->getCenterLocationOfBin(i) << "\t" << histogramNotDimer_[i] << "\n";
    }
    ofs.close();
}

void DistributionGivenDimer::printNumDimerPerResiduePerIter(std::ofstream& ofs)
{
    for (int i=0;i<DimerPerResidue_.size();i++)
    {
        ofs << DimerPerResidue_[i] << "\t";
    }
    
    ofs << "\n";
}

void DistributionGivenDimer::printNumDimerPerIter(std::ofstream& ofs)
{
    int time = simstate_.getTime();

    ofs << time << "\t" << numDimersPerIter_ << "\n";
}

void DistributionGivenDimer::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<histogram_.size();i++)
    {
        ofs << bin_ -> getCenterLocationOfBin(i) << "\t" << histogram_[i] << "\n";
    }

    ofs.close();
}
