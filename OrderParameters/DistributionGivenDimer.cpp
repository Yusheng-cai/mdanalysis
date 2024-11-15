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
    numatoms_ = res[0].atoms_.size();
    AngleWithSurface_.resize(numres_,0.0);
    COM_.resize(numres_);
    COMdistance_.resize(numres_);
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

    // read the COM indices for calculating distances between molecules
    DistanceCOM_.resize(numatoms_,0.0);
    std::iota(DistanceCOM_.begin(), DistanceCOM_.end(), 1);
    pack_.ReadVectorNumber("distanceCOM", ParameterPack::KeyType::Optional, DistanceCOM_);
    for (int i=0;i<DistanceCOM_.size();i++)
    {
        DistanceCOM_[i] -= 1;
    }

    // read in the cosine theta max and r max 
    pack_.ReadNumber("cosmax", ParameterPack::KeyType::Required, cosmax_);
    pack_.ReadNumber("rmax", ParameterPack::KeyType::Required, rmax_);

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});
    registerOutputFunction("histogramnotdimer", [this](std::string name) -> void {this->printHistogramNotdimer(name);});
    registerPerIterOutputFunction("dimers", [this](std::ofstream& ofs) -> void {this -> printNumDimerPerIter(ofs);});
    registerPerIterOutputFunction("dimerbetafactor", [this](std::ofstream& ofs) -> void {this -> printDimerBetaFactor(ofs);});
    registerPerIterOutputFunction("notdimerbetafactor", [this](std::ofstream& ofs) -> void {this -> printNotDimerBetaFactor(ofs);});
    registerOutputFileOutputs("dimers", [this](void) -> Real {return this->getNumDimers();});
    registerPerIterOutputFunction("dimerperres", [this](std::ofstream& ofs) -> void {this -> printNumDimerPerResiduePerIter(ofs);});
}

void DistributionGivenDimer::calculate()
{
    auto& res = getResidueGroup(resName_).getResidues();
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COMdistance_[i] = CalculationTools::getCOM(res[i], simstate_, DistanceCOM_);

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
    InsideIndices_ = InsidePVIndices(COM_, OutsideIndices_);
    int size = InsideIndices_.size();
    DimerPerResidue_.clear();
    DimerPerResidue_.resize(COM_.size(),0.0);
    numDimersPerIter_ = 0;

    #pragma omp parallel
    {
        std::vector<Real> dbuffer(COM_.size(),0.0);
        int numdimers = 0;
        #pragma omp for
        for (int i=0;i<size;i++)
        {
            for (int j=0;j<COM_.size();j++)
            {
                int index1 = InsideIndices_[i];

                if (index1 != j)
                {
                    Real3 distance;
                    Real distsq;
                    Real pairDistance;
                    Real costhetaPair;

                    simstate_.getSimulationBox().calculateDistance(COMdistance_[index1], COMdistance_[j], distance, distsq);

                    pairDistance = std::sqrt(distsq);
                    Real costheta = LinAlg3x3::DotProduct(uij_[index1], uij_[j]);
                    costhetaPair = costheta;

                    if (pairDistance <= rmax_ && costhetaPair <= cosmax_)
                    {
                        dbuffer[index1] += 1;
                        dbuffer[j]      += 1;
                        numdimers += 1;
                    }
                }
            }
        }

        #pragma omp critical
        {
            numDimersPerIter_ += numdimers;
        }

        #pragma omp critical
        {
            for (int i=0;i<COM_.size();i++)
            {
                DimerPerResidue_[i] += dbuffer[i];
            }
        }
    }

    // bin it into histograms
    for (int i=0;i<DimerPerResidue_.size();i++)
    {
        Real angle = AngleWithSurface_[i];

        if (bin_ -> isInRange(angle))
        {
            int binnum = bin_ -> findBin(angle);

            if (DimerPerResidue_[i] > 0)
            {
                histogram_[binnum] += 1;
            }
            else
            {
                histogramNotDimer_[binnum] += 1;
            }
        }
    }
}

void DistributionGivenDimer::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();
    Real sum1=0.0;
    Real sum2=0.0;
    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= numFrames;
        histogramNotDimer_[i] /= numFrames;

        sum1 += histogram_[i];
        sum2 += histogramNotDimer_[i];
    }

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= sum1;
        histogramNotDimer_[i] /= sum2;
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

void DistributionGivenDimer::printDimerBetaFactor(std::ofstream& ofs)
{
    auto& res = getResidueGroup(resName_).getResidues();
    int timestep = simstate_.getStep();
    ofs <<  timestep << " ";

    for (int i=0;i<DimerPerResidue_.size();i++)
    {
        if (DimerPerResidue_[i] > 0)
        {
            for (int j=0;j<res[i].atoms_.size();j++)
            {
                ofs << res[i].atoms_[j].atomNumber_ - 1 << " ";
            }
        }
    }

    ofs << "\n";
}

void DistributionGivenDimer::printNotDimerBetaFactor(std::ofstream& ofs)
{
    auto& res = getResidueGroup(resName_).getResidues();

    int timestep = simstate_.getStep(); 
    ofs << timestep << " ";

    for (int i=0;i<DimerPerResidue_.size();i++)
    {
        if (DimerPerResidue_[i] == 0)
        {
            for (int j=0;j<res[i].atoms_.size();j++)
            {
                ofs << res[i].atoms_[j].atomNumber_ -1 << " ";
            }
        }
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
