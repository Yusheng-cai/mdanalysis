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
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);

    bins_ = Binptr(new Bin(*binPack));
    numCountsPerBin_.resize(bins_->getNumbins());
    std::fill(numCountsPerBin_.begin(), numCountsPerBin_.end(),0);

    // add the residue groups
    addResidueGroup(resname1_);
    addResidueGroup(resname2_);

    auto& res1 = getResidueGroup(resname1_).getResidues();
    auto& res2 = getResidueGroup(resname2_).getResidues();
    COMIndices1_.resize(res1.size());
    COMIndices2_.resize(res2.size());
    std::iota(COMIndices1_.begin(), COMIndices1_.end(),0);
    std::iota(COMIndices2_.begin(), COMIndices2_.end(),0);

    input.pack_.ReadVectorNumber("comIndices1", ParameterPack::KeyType::Optional, COMIndices1_);
    input.pack_.ReadVectorNumber("comIndices2", ParameterPack::KeyType::Optional, COMIndices2_);

    ASSERT((res1.size() == res2.size()), "The size of residue 1 does not match that of residue 2");
}

void RDFresidue::calculate()
{
    const auto& res1 = getResidueGroup(resname1_).getResidues();
    const auto& res2 = getResidueGroup(resname2_).getResidues();

    distance_.clear();
    distanceBuffer_.clearBuffer();
    distanceBuffer_.set_master_object(distance_);

    COM1_.clear();
    COM1_.resize(res1.size());
    COM2_.clear();
    COM2_.resize(res2.size());

    // The size is half of N1*N2
    // distance_.resize(res1.size()*res2.size()*0.5);

    #pragma omp parallel for
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

    for (int i=0;i<distance_.size();i++)
    {
        int binNum = bins_->findBin(distance_[i]);
        numCountsPerBin_[binNum] += 1;
    }
}