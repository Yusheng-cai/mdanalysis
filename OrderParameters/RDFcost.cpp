#include "RDFcost.h"

RDFcost::RDFcost(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
    input.pack_.ReadNumber("headIndex", ParameterPack::KeyType::Required, headIndex_);
    input.pack_.ReadNumber("tailIndex", ParameterPack::KeyType::Required, tailIndex_);

    addResidueGroup(resName_);

    auto res = getResidueGroup(resName_).getResidues();

    COMIndices_.resize(res[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1);

    input.pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);

    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = Binptr(new Bin(*binPack));

    gcostheta_.resize(bin_->getNumbins());
    std::fill(gcostheta_.begin(), gcostheta_.end(), 0.0);

    // resize the atomic director vector
    uij_.resize(res.size());
}


void RDFcost::calculate()
{
    auto& res = getResidueGroup(resName_).getResidues();

    // first find the center of mass as well as the uij
    std::vector<Real3> COM(res.size());
    for (int i=0;i<res.size();i++)
    {
        auto& residue = res[i];

        Real3 com = CalculationTools::getCOM(residue, simstate_, COMIndices_);
        COM[i] = com;

        Real3 distance;
        Real dist_sq;

        Real3 positionHead = residue.atoms_[headIndex_].positions_;
        Real3 positionTail = residue.atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(positionHead, positionTail, distance, dist_sq);

        Real3 dist_norm = Qtensor::normalize_director(distance);

        uij_[i] = dist_norm;
    }

    // the number of pairs per bin in each iteration
    std::vector<int> NumPerBinIter_(bin_->getNumbins(),0);
    for (int i=0;i<res.size();i++)
    {
        for (int j=i+1;j<res.size();j++)
        {
            Real3 distance;
            Real dist_sq;

            simstate_.getSimulationBox().calculateDistance(COM[i], COM[j], distance, dist_sq);
            Real dist = std::sqrt(dist_sq);
            
            if (bin_ -> isInRange(dist))
            {
                int binNum = bin_ ->findBin(dist);
                Real dot_product = Qtensor::vec_dot(uij_[i], uij_[j]);

                gcostheta_[binNum] += dot_product;
            }
        }
    }
}