#include "Pcostz.h"

namespace CalculationRegistry
{
    static registry_<Pcostz> registerPcostz("pcostz");
}

Pcostz::Pcostz(const CalculationInput& input)
:Calculation(input)
{
    auto zbinPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Required);
    auto costBinPack = input.pack_.findParamPack("tbin", ParameterPack::KeyType::Required);
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,tailIndex_);
    input.pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    input.pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    bool outputRead = input.pack_.ReadString("output", ParameterPack::KeyType::Optional,outputName_);

    if(outputRead)
    {
        ofs_.open(outputName_);
    }

    headIndex_--;
    tailIndex_--;

    directionIndex_ = MapdirectionToIndex_.find(direction_) -> second;

    // add the residue group
    addResidueGroup(residueGroupName_);

    costBin_ = Binptr(new Bin(*costBinPack));
    zBin_    = Binptr(new Bin(*zbinPack));

    auto& res = getResidueGroup(residueGroupName_).getResidues();
    COMIndices_.resize(res[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1);
    input.pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);

    for (int i=0;i<COMIndices_.size();i++)
    {
        COMIndices_[i] -= 1;
    }

    COM_.resize(res.size());

    std::vector<Real> Zero(costBin_->getNumbins(),0);
    histogram2d_.resize(zBin_->getNumbins(), Zero);

    // resize the molecular director size
    uij_.resize(res.size());
}

void Pcostz::calculate()
{
    auto& res = getResidueGroup(residueGroupName_).getResidues();

    // find all the COM of the residues in the system
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

    // find the distances between all pairs of residues
    for (int i=0;i<res.size();i++)
    {
        Real dist = COM_[i][directionIndex_];
        Real cost = Qtensor::vec_dot(uij_[i], arr_);
        if (zBin_->isInRange(dist) && costBin_->isInRange(cost))
        {  
            int zbinNum = zBin_->findBin(dist);
            int tbinNum = costBin_->findBin(cost);


            histogram2d_[zbinNum][tbinNum] += 1;
        }
    }
}

void Pcostz::finishCalculate()
{
    int Numframes = simstate_.getTotalFrames();

    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            histogram2d_[i][j] /= Numframes;
        }
    }
}

void Pcostz::printOutput()
{
    if (ofs_.is_open())
    {
        ofs_ << std::fixed << std::setprecision(precision_);

        for (int i=0;i<histogram2d_.size();i++)
        {
            for (int j=0;j<histogram2d_[0].size();j++)
            {
                ofs_ << i << "\t" << j << "\t" << histogram2d_[i][j] << "\n";
            }
        }
        ofs_.close();
    }
}