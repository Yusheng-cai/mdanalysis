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
    input.pack_.ReadVectorNumber("comIndices", ParameterPack::KeyType::Optional, COMIndices_);

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

    // add the residue group to the system
    addResidueGroup(residueName_);

    // resize the binned matrix to number of bins
    BinnedMatrix_.resize(bin_->getNumbins());

    // resize the P2 to be number of bins size
    P2_.resize(bin_->getNumbins());
    std::fill(P2_.begin(), P2_.end(), 0.0);
}

void QtensorZ::calculate()
{
    // obtain the residue group by its name
    const auto& res = getResidueGroup(residueName_).getResidues();

    // obtain the COM
    std::vector<Real3> COM = CalculationTools::getCOM(res,simstate_);

    // make all the matrix in vector zero
    Matrix zeroMatrix;
    zeroMatrix.fill({});
    std::fill(BinnedMatrix_.begin(), BinnedMatrix_.end(), zeroMatrix);

    // for (int i=0;i<COM.size();i++)
    // {
    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << COM[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    std::vector<int> NumResPerBin_;
    NumResPerBin_.resize(bin_->getNumbins());
    std::fill(NumResPerBin_.begin(), NumResPerBin_.end(),0);

    // Bin the COMs
    for (int i=0;i<COM.size();i++)
    {
        Real num = COM[i][index_];

        int binNum = bin_ -> findBin(num);

        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

        Real3 diff;
        Real diff_sq;

        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, diff, diff_sq);

        // normalize the director and calculate the dyad
        Qtensor::normalize_director(diff);
        Matrix dyad = Qtensor::vec_dyadic(diff, diff);

        Qtensor::matrix_mult_inplace(dyad, 3.0);
        Matrix Qlocal = Qtensor::matrix_sub(dyad, Qtensor::matrix_Identity());

        Qtensor::matrix_accum_inplace(BinnedMatrix_[binNum], Qlocal);

        NumResPerBin_[binNum] += 1;
    }

    for (int i=0;i<BinnedMatrix_.size();i++)
    {
        Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 1.0/(2.0*NumResPerBin_[i]));

        auto result = Qtensor::OP_Qtensor(BinnedMatrix_[i]);

        P2_[i] += result.first;
    }
}

void QtensorZ::finishCalculate()
{
    int totalFrames = simstate_.getTotalFrames();

    for (int i=0;i<P2_.size();i++)
    {
        P2_[i] /= totalFrames;
    }
}