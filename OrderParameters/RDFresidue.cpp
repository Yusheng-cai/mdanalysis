#include "RDFresidue.h"

namespace CalculationRegistry
{
    registry_<RDFresidue> registerRDFresidue("RDFresidue");
}

RDFresidue::RDFresidue(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residue1", ParameterPack::KeyType::Required, resname1_);
    input.pack_.ReadString("residue2", ParameterPack::KeyType::Required, resname2_);
    bool readOutput = input.pack_.ReadString("output", ParameterPack::KeyType::Optional, outputName_);

    if (readOutput)
    {
        outputofs_.open(outputName_);
    }

    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);

    bins_ = Binptr(new Bin(*binPack));
    // initialize the volume 
    volume_.resize(bins_->getNumbins());
    Real dr = bins_ -> getStep();
    volume_[0] = 4.0/3.0*Constants::PI*std::pow(dr,3.0);

    for (int i=1;i<bins_->getNumbins();i++)
    {
        Real r = bins_ -> getLeftLocationOfBin(i);
        volume_[i] = 4.0*Constants::PI*std::pow(r,2.0)*dr;
    }

    // add the residue groups
    addResidueGroup(resname1_);
    addResidueGroup(resname2_);

    auto& res1 = getResidueGroup(resname1_).getResidues();
    auto& res2 = getResidueGroup(resname2_).getResidues();
    COMIndices1_.resize(res1[0].atoms_.size());
    COMIndices2_.resize(res2[0].atoms_.size());
    std::iota(COMIndices1_.begin(), COMIndices1_.end(),1);
    std::iota(COMIndices2_.begin(), COMIndices2_.end(),1);
    
    input.pack_.ReadVectorNumber("comIndices1", ParameterPack::KeyType::Optional, COMIndices1_);
    input.pack_.ReadVectorNumber("comIndices2", ParameterPack::KeyType::Optional, COMIndices2_);

    ASSERT((res1.size() == res2.size()), "The size of residue 1 does not match that of residue 2");

    // user input should be 1 based
    for (int i=0;i<COMIndices1_.size();i++)
    {
        COMIndices1_[i] -=1;
        COMIndices2_[i] -=1;
    }


    // resize rdf vector to be same lengths as number of bins
    rdf_.resize(bins_->getNumbins());
}

void RDFresidue::calculate()
{
    const auto& res1 = getResidueGroup(resname1_).getResidues();
    const auto& res2 = getResidueGroup(resname2_).getResidues();

    distance_.clear();
    distanceBuffer_.clearBuffer();
    distanceBuffer_.set_master_object(distance_);

    int N = res1.size();

    // obtain the density (num Residues/ nm3)
    Real rho = res1.size()/simstate_.getSimulationBox().getVolume();

    COM1_.clear();
    COM1_.resize(res1.size());
    COM2_.clear();
    COM2_.resize(res2.size());


    for (int i=0;i<COM1_.size();i++)
    {
        Real3 C1 = CalculationTools::getCOM(res1[i], simstate_, COMIndices1_);
        Real3 C2 = CalculationTools::getCOM(res2[i], simstate_, COMIndices2_);

        COM1_[i] = C1;
        COM2_[i] = C2;
    }


    #pragma omp parallel
    {
        auto& buffer = distanceBuffer_.access_buffer_by_id();
        #pragma omp for
        for (int i=0;i<COM1_.size();i++)
        {
            for (int j=i+1;j<COM2_.size();j++)
            {
                Real3 distance;
                Real sq_dist;
                simstate_.getSimulationBox().calculateDistance(COM1_[i], COM2_[j],distance, sq_dist);

                buffer.push_back(std::sqrt(sq_dist));
            }
        }
    }


    int size = distance_.size();
    for (auto it = distanceBuffer_.beginworker();it != distanceBuffer_.endworker(); it++)
    {
        size += it -> size();
    }

    distance_.reserve(size);

    for (auto it = distanceBuffer_.beginworker();it != distanceBuffer_.endworker();it++)
    {
        distance_.insert(distance_.end(), it -> begin(), it -> end());
    }

    std::vector<int> NumCountsPerBin(bins_->getNumbins(), 0);

    #pragma omp parallel
    {    
        std::vector<int> NumberCountsLocal(bins_->getNumbins(),0);

        #pragma omp for
        for (int i=0;i<distance_.size();i++)
        {
            if (bins_->isInRange(distance_[i]))
            {
                int binNum = bins_->findBin(distance_[i]);

                NumberCountsLocal[binNum] += 1;
            }
        }

        #pragma omp critical
        {
            for (int i=0;i<bins_->getNumbins();i++)
            {
                NumCountsPerBin[i] += NumberCountsLocal[i];
            }
        }
    }

    for (int i=0;i<rdf_.size();i++)
    {
        rdf_[i] = rdf_[i] + 2.0*NumCountsPerBin[i]/(rho*volume_[i]*N); 
    }

}

void RDFresidue::finishCalculate()
{
    // The differential r is the step in the bins
    Real dr = bins_->getStep();
    int numBins = bins_->getNumbins();

    // normalize the numCountsPerBin By the number of frames first
    int numFrames_ = simstate_.getTotalFrames();

    for (int i=0;i<rdf_.size();i++)
    {
        rdf_[i] = rdf_[i]/numFrames_;
    }
}

void RDFresidue::printOutput()
{
    if (outputofs_.is_open())
    {
        outputofs_ << std::fixed << std::setprecision(precision_);
        outputofs_ << "# bins\trdf\n";

        for (int i=0;i<rdf_.size();i++)
        {
            outputofs_ << bins_->getLeftLocationOfBin(i) << " "; 
            outputofs_ << rdf_[i] << "\n";
        }
        outputofs_.close();
    }

}