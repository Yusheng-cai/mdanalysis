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
    COMIndex_provided_ = input.pack_.ReadVectorNumber("COMindices", ParameterPack::KeyType::Optional, COMIndices_);

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
    COM_.clear();
    COM_.resize(res.size());

    for (int i=0;i<res.size();i++)
    {
        int residueSize = res[i].atoms_.size();
        Real3 COMperAtom_;
        if (! COMIndex_provided_)
        {
            std::vector<int> chosen_indices_;
            chosen_indices_.resize(residueSize);
            std::iota(chosen_indices_.begin(), chosen_indices_.end(), 0);

            COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, chosen_indices_);
        }
        else
        {
            COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        }
        COM_[i] = COMperAtom_;

        std::cout << "COM for residue " << i << std::endl;
        for (int j=0;j<3;j++)
        {
            std::cout << COMperAtom_[j] << " ";
        }
        std::cout << "\n";
    }

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
    for (int i=0;i<COM_.size();i++)
    {
        Real num = COM_[i][index_];

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