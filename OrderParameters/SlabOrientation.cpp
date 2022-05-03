#include "SlabOrientation.h"

namespace CalculationRegistry
{
    static registry_<SlabOrientation> registerSlabOrientation("SlabOrientation");
}

SlabOrientation::SlabOrientation(const CalculationInput& input)
:Calculation(input)
{
    // initialize the bins first --> zbin and thetabin
    auto zbinPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);
    auto costBinPack = input.pack_.findParamPack("tbin", ParameterPack::KeyType::Required);
    costBin_ = Binptr(new Bin(*costBinPack));
    if(zbinPack != nullptr)
    {
        zBin_    = Binptr(new Bin(*zbinPack));
        Znumbins_= zBin_->getNumbins();
    }
    else
    {
        pack_.ReadNumber("Znumbins", ParameterPack::KeyType::Required, Znumbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Required, above_);
        zBin_    = Binptr(new Bin());
        usingMinMax_ = true;
    }

    // register the output functions
    RegisterOutputs();

    // read the inputs 
    ReadInputs();

    // initialize residue group
    initializeResidueGroup(residueGroupName_);
    auto& res = getResidueGroup(residueGroupName_);

    // size histogram 2d      
    histogram2d_.resize(Znumbins_, std::vector<Real>(costBin_->getNumbins(),0.0));

    // resize the molecular director size
    uij_.resize(res.size());

    // resize the number of residues to size of residue
    numResiduePerBin_.resize(res.size());
    std::fill(numResiduePerBin_.begin(), numResiduePerBin_.end(),0);
}

void SlabOrientation::RegisterOutputs()
{
    registerOutputFunction("PCosthetaZ", [this](std::string name)->void {this -> printHistogram(name);});
}

void SlabOrientation::ReadInputs()
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,tailIndex_);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    // correct for the indexing on head and tail index 
    headIndex_--;
    tailIndex_--;

    // find the direction
    directionIndex_ = MapdirectionToIndex_.find(direction_) -> second;
}

void SlabOrientation::calculate()
{
    // obtain the residue group
    auto& res = getResidueGroup(residueGroupName_).getResidues();

    // find all the COM of the residues in the system
    for (int i=0;i<res.size();i++)
    {
        Real3 com = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = com;

        Real3 distance;
        Real dist_sq;
        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, dist_sq);

        Real3 normalized_dir = Qtensor::normalize_director(distance);
        uij_[i] = normalized_dir;
    }

    // bin using min max of the molecules if needed 
    if (usingMinMax_)
    {
        binUsingMinMax();
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
            numResiduePerBin_[zbinNum] += 1;
        }
    }
}

void SlabOrientation::binUsingMinMax()
{
    Real slight_shift=1e-3;
    std::vector<Real> ZdirectionNum;

    for (int i=0;i<COM_.size();i++)
    {
        Real val = COM_[i][directionIndex_];

        if (val > above_)
        {
            ZdirectionNum.push_back(val);
        }
    }

    auto maxit = std::max_element(ZdirectionNum.begin(), ZdirectionNum.end());
    auto minit = std::min_element(ZdirectionNum.begin(), ZdirectionNum.end());

    Real max = *maxit + slight_shift;
    Real min = *minit - slight_shift;
    Range2 range = {{min, max}};

    std::cout << "Max = " << max << " Min = " << min << std::endl;

    zBin_ -> update(range, Znumbins_);
}

void SlabOrientation::finishCalculate()
{
    int Numframes = simstate_.getTotalFrames();

    // first normalize by the number of frames 
    for (int i=0;i<numResiduePerBin_.size();i++)
    {
        numResiduePerBin_[i] /= Numframes;
    }

    // normalize the entire histogram such as \sum \sum P(z, cos(\beta)) = 1
    // first divide P(z, cos(\beta)) by numframes
    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            histogram2d_[i][j] /= Numframes;
        }
    }

    // Real sum = 0.0;
    // for (int i=0;i<histogram2d_.size();i++)
    // {
    //     for (int j=0;j<histogram2d_[0].size();j++)
    //     {
    //         sum += histogram2d_[i][j];
    //     }
    // }

    // for (int i=0;i<histogram2d_.size();i++)
    // {
    //     for (int j=0;j<histogram2d_[0].size();j++)
    //     {
    //         histogram2d_[i][j] = histogram2d_[i][j] / sum;
    //     }
    // }
}

void SlabOrientation::printHistogram(std::string name)
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
