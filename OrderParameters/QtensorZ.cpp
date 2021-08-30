#include "QtensorZ.h"

namespace CalculationRegistry
{
    registry_<QtensorZ> registerQtensorZ("qtensorZ");
}

QtensorZ::QtensorZ(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residuegroup", ParameterPack::KeyType::Required, residueName_);
    input.pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required, headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required, tailIndex_);
    input.pack_.ReadNumber("ignorelessthan", ParameterPack::KeyType::Optional, ignoreP2LessThan_);

    bool readP2z =input.pack_.ReadString("p2zOutput", ParameterPack::KeyType::Optional, p2ZOutput_);

    if (readP2z)
    {
        p2zofs_.open(p2ZOutput_);
    }

    // add the residue group to the system
    addResidueGroup(residueName_);

    auto& res = getResidueGroup(residueName_);

    // COM Indices start with 1
    COMIndices_.resize(res.getResidues()[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1.0);
    COMIndex_provided_ = input.pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);

    headIndex_ --;
    tailIndex_ --;

    for (int i=0;i<COMIndices_.size();i++)
    {
        COMIndices_[i] -= 1;
    }
 
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);

    // create the bin object
    bin_ = Binptr(new Bin(*binPack));

    auto it = MapNameToDirection.find(direction_);
    ASSERT((it != MapNameToDirection.end()), "The direction " << direction_ << " is not recognized.");
    index_ = it -> second;

    
    // resize the binned matrix to number of bins
    BinnedMatrix_.resize(bin_->getNumbins());

    // resize the P2 to be number of bins size
    P2_.resize(bin_->getNumbins());
    std::fill(P2_.begin(), P2_.end(), 0.0);

    // resize number of residues per bin
    NumResPerBin_.resize(bin_ -> getNumbins());
    std::fill(NumResPerBin_.begin(), NumResPerBin_.end(), 0);
}

void QtensorZ::calculate()
{
    // obtain the residue group by its name
    const auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());

    for (int i=0;i<res.size();i++)
    {
        int residueSize = res[i].atoms_.size();
        Real3 COMperAtom_;

        COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = COMperAtom_;
    }

    // make all the matrix in vector zero
    Matrix zeroMatrix;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            zeroMatrix[i][j] = 0;
        }
    }
    std::fill(BinnedMatrix_.begin(), BinnedMatrix_.end(), zeroMatrix);

    std::vector<int> NumResPerBinIter_;
    NumResPerBinIter_.resize(bin_->getNumbins());
    std::fill(NumResPerBinIter_.begin(), NumResPerBinIter_.end(),0);

    // Bin the COMs
    for (int i=0;i<COM_.size();i++)
    {
        Real num = COM_[i][index_];

        // only perform these operations of num is in range of bin
        if (bin_->isInRange(num))
        {
            int binNum = bin_ -> findBin(num);

            Real3 headPos = res[i].atoms_[headIndex_].positions_;
            Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

            Real3 diff;
            Real diff_sq;

            simstate_.getSimulationBox().calculateDistance(headPos, tailPos, diff, diff_sq);

            // normalize the director and calculate the dyad
            Real3 diffnorm = Qtensor::normalize_director(diff);
            Matrix dyad = Qtensor::vec_dyadic(diffnorm, diffnorm);

            Qtensor::matrix_mult_inplace(dyad, 3.0);
            Matrix Qlocal = Qtensor::matrix_sub(dyad, Qtensor::matrix_Identity());

            Qtensor::matrix_accum_inplace(BinnedMatrix_[binNum], Qlocal);

            NumResPerBinIter_[binNum] += 1;
            NumResPerBin_[binNum] += 1;
        }
    }

    int sum_ = 0;
    for (int i=0;i<NumResPerBinIter_.size();i++)
    {
        sum_ += NumResPerBinIter_[i];
    }

    for (int i=0;i<BinnedMatrix_.size();i++)
    {
        if (NumResPerBinIter_[i] != 0)
        {
            Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 1.0/(2.0*NumResPerBinIter_[i]));

            auto result = Qtensor::OP_Qtensor(BinnedMatrix_[i]);
            P2_[i] += result.first;
            std::cout << "Bin " << i << ", P2 = " << result.first << std::endl;
        }
    }
}

void QtensorZ::finishCalculate()
{
    int totalFrames = simstate_.getTotalFrames();

    for (int i=0;i<P2_.size();i++)
    {
        P2_[i] /= totalFrames;
        NumResPerBin_[i] /= totalFrames;

        if ( NumResPerBin_[i] < ignoreP2LessThan_ )
        {
            P2_[i] = 0.0;
        }
    }
}

void QtensorZ::printOutput()
{
    if (p2zofs_.is_open())
    {
        p2zofs_ << "#z\tp2\tNumber\n";

        for (int i=0;i<P2_.size();i++)
        {
            p2zofs_ << bin_->getLeftLocationOfBin(i) << "\t" << P2_[i] << "\t" << NumResPerBin_[i] << std::endl;
        }
        p2zofs_.close();
    }
}