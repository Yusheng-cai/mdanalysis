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
    input.pack_.ReadString("residuegroup", ParameterPack::KeyType::Required, residueGroupName_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,tailIndex_);
    input.pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    input.pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    input.pack_.ReadNumber("ignorelessthan", ParameterPack::KeyType::Optional, ignoreBelow_);

    registerOutputFunction("pcostz", [this](std::string name)->void {this -> printHistogram(name);});
    registerPerIterOutputFunction("pcostz", [this](std::ofstream& ofs)->void{this -> printHistogramPerIter(ofs);});

    headIndex_--;
    tailIndex_--;

    // find the direction
    directionIndex_ = MapdirectionToIndex_.find(direction_) -> second;
    costBin_ = Binptr(new Bin(*costBinPack));
    zBin_    = Binptr(new Bin(*zbinPack));

    // initialize residue group
    initializeResidueGroup(residueGroupName_);
    auto& res = getResidueGroup(residueGroupName_);

    std::vector<Real> Zero(costBin_->getNumbins(),0);
    histogram2d_.resize(zBin_->getNumbins(), Zero);

    // resize the molecular director size
    uij_.resize(res.size());

    // resize the number of residues to size of residue
    numResiduePerBin_.resize(res.size());
    std::fill(numResiduePerBin_.begin(), numResiduePerBin_.end(),0);
}

void Pcostz::calculate()
{
    auto& res = getResidueGroup(residueGroupName_).getResidues();

    numResiduePerBinIter_.clear();
    numResiduePerBinIter_.resize(numResiduePerBin_.size(), 0);

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

    histogramIter_.clear();
    histogramIter_.resize(histogram2d_.size(), std::vector<Real>(histogram2d_[0].size(), 0.0));

    // find the distances between all pairs of residues
    for (int i=0;i<res.size();i++)
    {
        Real dist = COM_[i][directionIndex_];
        Real cost = Qtensor::vec_dot(uij_[i], arr_);

        ASSERT((cost >= -1 && cost <= 1), "cosine(theta) is not within range of -1 and 1");
        if (zBin_->isInRange(dist) && costBin_->isInRange(cost))
        {  
            int zbinNum = zBin_->findBin(dist);
            int tbinNum = costBin_->findBin(cost);


            histogram2d_[zbinNum][tbinNum] += 1;
            histogramIter_[zbinNum][tbinNum] += 1;

            numResiduePerBin_[zbinNum] += 1;
            numResiduePerBinIter_[zbinNum] += 1;
        }
    }
}

void Pcostz::finishCalculate()
{
    int Numframes = simstate_.getTotalFrames();

    for (int i=0;i<numResiduePerBin_.size();i++)
    {
        numResiduePerBin_[i] /= Numframes;
    }

    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            histogram2d_[i][j] /= Numframes;

            if (numResiduePerBin_[i] < ignoreBelow_)
            {
                histogram2d_[i][j] = 0;
            }
        }
    }

    std::vector<Real> Rowsum_(histogram2d_.size(),0);

    for (int i=0;i<histogram2d_.size();i++)
    {
        Real sum = 0;
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            sum += histogram2d_[i][j];
        }
        Rowsum_[i] = sum;
    }

    // normalize each row by the sum if sum is not 0
    for (int i=0;i<histogram2d_.size();i++)
    {
        if (Rowsum_[i] != 0)
        {
            for (int j=0;j<histogram2d_[0].size();j++)
            {
                histogram2d_[i][j] /= Rowsum_[i];
            }
        }
    }
}

void Pcostz::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            ofs << i << "\t" << j << "\t" << histogram2d_[i][j] << "\n";
        }
    }
    ofs.close();
}

void Pcostz::printHistogramPerIter(std::ofstream& ofs)
{
    // if the file is open, then we perform the calculations
    int stepnum = simstate_.getFrameNumber();

    std::cout << "numResiduePerBinIter[1] = " << numResiduePerBinIter_[1] << std::endl; 
    for (int i=0;i<histogramIter_.size();i++)
    {
        for (int j=0;j<histogramIter_[0].size();j++)
        {
            if (numResiduePerBinIter_[i] < ignoreBelow_)
            {
                histogramIter_[i][j] = 0;
            }
        }
    }

    // perform row sum
    std::vector<Real> RowSum(histogramIter_.size(),0.0);
    for (int i=0;i<histogramIter_.size();i++)
    {
        Real sum_ = 0.0;
        for (int j=0;j<histogramIter_[0].size();j++)
        {
            sum_ += histogramIter_[i][j];
        }
        RowSum[i] = sum_;
    }
    
    // normalize each row
    for (int i=0;i<histogramIter_.size();i++)
    {
        if (RowSum[i] > 0)
        {
            for (int j=0;j<histogramIter_[0].size();j++)
            {
                histogramIter_[i][j] /= RowSum[i];
            }
        } 
    }

    for (int i=0;i<histogramIter_.size();i++)
    {
        for (int j=0;j<histogramIter_[0].size();j++)
        {
            ofs << histogramIter_[i][j] << " ";
        }
    }

    ofs << "\n";
}
