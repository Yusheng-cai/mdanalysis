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

    registerOutputFunction("p2z", [this](std::string name) -> void{this -> printP2z(name);});
    registerPerIterOutputFunction("p2z", [this](std::ofstream& ofs) -> void {this -> printPerIterP2z(ofs);});
    registerPerIterOutputFunction("Qtensor", [this](std::ofstream& ofs) -> void {this -> printPerIterQtensor(ofs);});
    registerPerIterOutputFunction("ev", [this](std::ofstream& ofs) -> void{this -> printPerIterev(ofs);});

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
    // make all the matrix in vector zero
    Matrix zeroMatrix = {};
    std::fill(BinnedMatrix_.begin(), BinnedMatrix_.end(), zeroMatrix);
    BinnedMatrix_.resize(bin_->getNumbins(),zeroMatrix);

    // resize the P2 to be number of bins size
    P2_.resize(bin_->getNumbins());
    std::fill(P2_.begin(), P2_.end(), 0.0);

    // resize number of residues per bin
    NumResPerBin_.resize(bin_ -> getNumbins());
    std::fill(NumResPerBin_.begin(), NumResPerBin_.end(), 0);

    // eigenvector
    std::array<Real,3> zeroArray = {{0,0,0}};
    eigvec_.resize( bin_ -> getNumbins());
    std::fill(eigvec_.begin(), eigvec_.end(), zeroArray);
}

void QtensorZ::calculate()
{
    // obtain the residue group by its name
    const auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());

    // initialize Per iter items
    evPerIter_.clear();
    P2PerIter_.clear();
    BinnedMatrixIter_.clear();
    NumResPerBinIter_.clear();

    evPerIter_.resize(bin_ -> getNumbins(), {{0,0,0}});
    P2PerIter_.resize(bin_ -> getNumbins(), 0.0);
    // make all the matrix in vector zero
    Matrix zeroMatrix;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            zeroMatrix[i][j] = 0;
        }
    }
    BinnedMatrixIter_.resize(bin_ -> getNumbins(), zeroMatrix);
    NumResPerBinIter_.resize(bin_ ->getNumbins(), 0.0);

    // obtain the center of mass of each of the residues
    for (int i=0;i<res.size();i++)
    {
        int residueSize = res[i].atoms_.size();
        Real3 COMperAtom_;

        COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = COMperAtom_;
    }

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

            Qtensor::matrix_accum_inplace(BinnedMatrixIter_[binNum], Qlocal);

            NumResPerBinIter_[binNum] += 1;
            NumResPerBin_[binNum] += 1;
        }
    }


    for (int i=0;i<BinnedMatrixIter_.size();i++)
    {
        if (NumResPerBinIter_[i] != 0)
        {
            Qtensor::matrix_mult_inplace(BinnedMatrixIter_[i], 1.0/(2.0*NumResPerBinIter_[i]));

            auto result = Qtensor::OP_Qtensor(BinnedMatrixIter_[i]);
            auto ev     = Qtensor::orderedeig_Qtensor(BinnedMatrixIter_[i]).first;

            P2_[i] += result.first;
            P2PerIter_[i] = result.first;

            for (int k=0;k<3;k++)
            {
                Real value = std::pow(ev[k][0],2.0);
                eigvec_[i][k] += value;
                evPerIter_[i][k] = value;
            }
        }
    }

    for (int i=0;i<BinnedMatrixIter_.size();i++)
    {
        Qtensor::matrix_accum_inplace(BinnedMatrix_[i], BinnedMatrixIter_[i]);
    }
}

void QtensorZ::printOutputOnStep()
{
    // first make those p2 which has less residue than specified to be 0
    for (int i=0;i<NumResPerBinIter_.size();i++)
    {
        if (NumResPerBinIter_[i] < ignoreP2LessThan_)
        {
            P2PerIter_[i] = 0;
            evPerIter_[i] = {{0,0,0}};
            NumResPerBinIter_[i] = 0;
            BinnedMatrixIter_[i] = zeroMatrix_;
        }
    }
}

void QtensorZ::printP2z(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    ofs << "#z\tp2\tNumber\tQxx\tQxy\tQxz\tQyy\tQyz\tQzz\tnx\tny\tnz\tp2avg\n";

    for (int i=0;i<P2_.size();i++)
    {
        ofs << bin_->getLeftLocationOfBin(i) << "\t" << P2_[i] << "\t" << NumResPerBin_[i];
        ofs << "\t" << BinnedMatrix_[i][0][0] << "\t" << BinnedMatrix_[i][0][1] << "\t" << BinnedMatrix_[i][0][2]; 
        ofs << "\t" << BinnedMatrix_[i][1][1] << "\t" << BinnedMatrix_[i][1][2];
        ofs << "\t" << BinnedMatrix_[i][2][2];
        ofs << "\t" << eigvec_[i][0] << "\t" << eigvec_[i][1] << "\t" << eigvec_[i][2] << "\t";
        ofs << P2avg_[i] << "\n";
    }
    ofs.close();
}

void QtensorZ::printPerIterP2z(std::ofstream& ofs)
{
    for (int i=0;i<P2PerIter_.size();i++)
    {
        ofs << P2PerIter_[i] << " ";
    } 
    ofs << "\n";
}

void QtensorZ::printPerIterNum(std::ofstream& ofs)
{
    for (int i=0;i<NumResPerBinIter_.size();i++)
    {
        ofs << NumResPerBinIter_[i] << " ";
    }
    ofs << "\n";
}

void QtensorZ::printPerIterQtensor(std::ofstream& ofs)
{
    for (int i=0;i<BinnedMatrixIter_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                perIterQtensorofs_ << BinnedMatrixIter_[i][j][k] << " ";
            }
        }
    }
    perIterQtensorofs_ << "\n";
}

void QtensorZ::printPerIterev(std::ofstream& ofs)
{
    for (int i=0;i<evPerIter_.size();i++)
    {
        ofs << evPerIter_[i][0] << " " << evPerIter_[i][1] << " " << evPerIter_[i][2] << " ";
    }
    ofs << "\n";
}

void QtensorZ::finishCalculate()
{
    int totalFrames = simstate_.getTotalFrames();
    std::cout << "Total frame = " << totalFrames << std::endl;

    // average the p2 over time
    for (int i=0;i<P2_.size();i++)
    {
        P2_[i] /= totalFrames;
        NumResPerBin_[i] /= totalFrames;

        if ( NumResPerBin_[i] < ignoreP2LessThan_ )
        {
            P2_[i] = 0.0;
        }
    }

    // average the Qtensor over time
    P2avg_.resize(bin_->getNumbins());
    for (int i=0;i<BinnedMatrix_.size();i++)
    {
        Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 1.0/totalFrames);
        if (NumResPerBin_[i] < ignoreP2LessThan_)
        {
            Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 0.0);
        }
        auto ans = Qtensor::OP_Qtensor(BinnedMatrix_[i]);
        P2avg_[i] = ans.first;
    }

    // average the eigenvectors over time
    for (int i=0;i<eigvec_.size();i++)
    {
        auto ev = Qtensor::orderedeig_Qtensor(BinnedMatrix_[i]).first;
        for (int j=0;j<3;j++)
        {
            eigvec_[i][j] = ev[j][0];
        }

        if (NumResPerBin_[i] == 0 || NumResPerBin_[i] < ignoreP2LessThan_)
        {
            eigvec_[i] = {{0,0,0}};
        }
    }

    // finally, make the number of residues per bin to be 0
    for (int i=0;i<NumResPerBin_.size();i++)
    {
        if (NumResPerBin_[i] < ignoreP2LessThan_)
        {
            NumResPerBin_[i] = 0;
        }
    }
}