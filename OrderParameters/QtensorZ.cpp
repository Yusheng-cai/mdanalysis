#include "QtensorZ.h"

namespace CalculationRegistry
{
    registry_<QtensorZ> registerQtensorZ("qtensorZ");
}

QtensorZ::QtensorZ(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    input.pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailIndex_);
    input.pack_.ReadNumber("ignorelessthan", ParameterPack::KeyType::Optional, ignoreP2LessThan_);

    registerOutputFunction("p2z", [this](std::string name) -> void{this -> printP2z(name);});
    registerPerIterOutputFunction("p2z", [this](std::ofstream& ofs) -> void {this -> printPerIterP2z(ofs);});
    registerPerIterOutputFunction("Qtensor", [this](std::ofstream& ofs) -> void {this -> printPerIterQtensor(ofs);});
    registerPerIterOutputFunction("ev", [this](std::ofstream& ofs) -> void{this -> printPerIterev(ofs);});
    registerPerIterOutputFunction("num", [this](std::ofstream& ofs) -> void{this -> printPerIterNum(ofs);});
    registerPerIterOutputFunction("p2zbeta", [this](std::ofstream& ofs) -> void {this ->printPerIterP2zBeta(ofs);});
    registerPerIterOutputFunction("evbeta", [this](std::ofstream& ofs) -> void {this -> printPerItereveczBeta(ofs);});

    // add the residue group to the system
    initializeResidueGroup(residueName_);
    headIndex_--;
    tailIndex_--;
    
    // initialize bin
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Optional);
    // create the bin object , if specified, then used user specified bin, else we use min max of the COM
    if (binPack != nullptr)
    {
        bin_ = Binptr(new Bin(*binPack));
        numbins_ = bin_ -> getNumbins();
    }
    else
    {
        input.pack_.ReadNumber("numbins", ParameterPack::KeyType::Required, numbins_);
        input.pack_.ReadNumber("above", ParameterPack::KeyType::Optional, above_);
        bin_ = Binptr(new Bin());
        usingMinMaxBins_ = true;
    }

    // initialize direction
    auto it = MapNameToDirection.find(direction_);
    ASSERT((it != MapNameToDirection.end()), "The direction " << direction_ << " is not recognized.");
    index_ = it -> second;
    
    // resize the binned matrix to number of bins
    // make all the matrix in vector zero
    Matrix zeroMatrix = {};
    std::fill(BinnedMatrix_.begin(), BinnedMatrix_.end(), zeroMatrix);
    BinnedMatrix_.resize(numbins_,zeroMatrix);

    // resize the P2 to be number of bins size
    P2_.resize(numbins_);
    std::fill(P2_.begin(), P2_.end(), 0.0);

    // resize number of residues per bin
    NumResPerBin_.resize(numbins_);
    std::fill(NumResPerBin_.begin(), NumResPerBin_.end(), 0);

    // eigenvector
    std::array<Real,3> zeroArray = {{0,0,0}};
    eigvec_.resize( numbins_);
    std::fill(eigvec_.begin(), eigvec_.end(), zeroArray);

    auto& res = getResidueGroup(residueName_);

    // res index to bin index 
    ResIndexToBinIndex_.resize(res.size(),0);
    BetaFactors_.resize(res.getAtomSize(),0.0);
}

void QtensorZ::binUsingMinMax()
{
    std::vector<Real> zdir;
    for (int i=0;i<COM_.size();i++)
    {
        if (COM_[i][index_] > above_)
        {
            zdir.push_back(COM_[i][index_]);
        }
    }

    auto minit = std::min_element(zdir.begin(), zdir.end());
    auto maxit = std::max_element(zdir.begin(), zdir.end());

    Real min   = *minit;
    Real max   = *maxit;

    Range range = {{min, max}};

    bin_->update(range, numbins_);
    std::cout << "Range = {" <<  min << " " << max << "}." << std::endl;
    std::cout << "step size = " << bin_ -> getStep() << std::endl;
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
    evPerIter_.resize(numbins_, {{0,0,0}});
    P2PerIter_.resize(numbins_, 0.0);
    ResIndexToBinIndex_.resize(res.size(),0);

    // make all the matrix in vector zero
    Matrix zeroMatrix;
    zeroMatrix.fill({});
    BinnedMatrixIter_.resize(numbins_, zeroMatrix);
    NumResPerBinIter_.resize(numbins_, 0.0);

    // obtain the center of mass of each of the residues
    for (int i=0;i<res.size();i++)
    {
        Real3 COMperAtom_;
        COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = COMperAtom_;
    }

    // update the bins using COM min and max 
    if (usingMinMaxBins_)
    {
        binUsingMinMax();
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

            ResIndexToBinIndex_[i] = binNum;

            Real3 diff;
            Real diff_sq;

            simstate_.getSimulationBox().calculateDistance(headPos, tailPos, diff, diff_sq);

            // normalize the director and calculate the dyad
            Real3 diffnorm = Qtensor::normalize_director(diff);
            Matrix Qlocal  = LinAlg3x3::LocalQtensor(diffnorm);

            Qtensor::matrix_accum_inplace(BinnedMatrixIter_[binNum], Qlocal);

            NumResPerBinIter_[binNum] += 1;
            NumResPerBin_[binNum] += 1;
        }
    }

    // obtain the per iter values
    for (int i=0;i<BinnedMatrixIter_.size();i++)
    {
        if (NumResPerBinIter_[i] != 0)
        {
            Qtensor::matrix_mult_inplace(BinnedMatrixIter_[i], 1.0/(2.0*NumResPerBinIter_[i]));

            auto result = Qtensor::OP_Qtensor(BinnedMatrixIter_[i]);
            auto ev     = Qtensor::orderedeig_Qtensor(BinnedMatrixIter_[i]).first;

            P2PerIter_[i] = result.first;

            for (int k=0;k<3;k++)
            {
                Real value = ev[k][0];
                evPerIter_[i][k] = value;
            }
        }
    }

    // sum up the matrices from per iter to the total one
    for (int i=0;i<BinnedMatrixIter_.size();i++)
    {
        Qtensor::matrix_accum_inplace(BinnedMatrix_[i], BinnedMatrixIter_[i]);
    }

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

void QtensorZ::printPerIterP2zBeta(std::ofstream& ofs)
{
    const auto& res = getResidueGroup(residueName_).getResidues();

    BetaFactors_.resize(getResidueGroup(residueName_).getAtomSize(),0.0);
    int index=0;

    // Now we calculate the beta factors 
    for (int i=0;i<res.size();i++)
    {
        int binNum = ResIndexToBinIndex_[i];
        Real binP2 = P2PerIter_[binNum];

        for (int j=0;j<res[i].atoms_.size();j++)
        {
            BetaFactors_[index] = binP2;
            index ++;
        }
    }

    int framenum = simstate_.getFrameNumber();

    ofs << framenum << " ";

    for (int i=0;i<BetaFactors_.size();i++)
    {
        ofs << BetaFactors_[i] << " ";
    }
    ofs << "\n";
}

void QtensorZ::printPerItereveczBeta(std::ofstream& ofs)
{
    const auto& res = getResidueGroup(residueName_).getResidues();

    std::vector<Real> values(getResidueGroup(residueName_).getAtomSize(),0.0);
    int index=0;

    // Now we calculate the beta factors 
    for (int i=0;i<res.size();i++)
    {
        int binNum = ResIndexToBinIndex_[i];
        Real3 ev = evPerIter_[binNum];
        Real  val= std::pow(ev[2],2); 

        for (int j=0;j<res[i].atoms_.size();j++)
        {
            values[index] = val;
            index ++;
        }
    }

    int framenum = simstate_.getFrameNumber();

    ofs << framenum << " ";

    for (int i=0;i<values.size();i++)
    {
        ofs << values[i] << " ";
    }
    ofs << "\n";
}


void QtensorZ::printP2z(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    ofs << "#z\tp2\tNumber\tQxx\tQxy\tQxz\tQyy\tQyz\tQzz\tnx\tny\tnz\n";

    for (int i=0;i<P2avg_.size();i++)
    {
        ofs << bin_->getLeftLocationOfBin(i) << "\t" << P2avg_[i] << "\t" << NumResPerBin_[i];
        ofs << "\t" << BinnedMatrix_[i][0][0] << "\t" << BinnedMatrix_[i][0][1] << "\t" << BinnedMatrix_[i][0][2]; 
        ofs << "\t" << BinnedMatrix_[i][1][1] << "\t" << BinnedMatrix_[i][1][2];
        ofs << "\t" << BinnedMatrix_[i][2][2];
        ofs << "\t" << eigvec_[i][0] << "\t" << eigvec_[i][1] << "\t" << eigvec_[i][2] << "\n";
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
                ofs << BinnedMatrixIter_[i][j][k] << " ";
            }
        }
    }
    ofs << "\n";
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

    // average the Qtensor over time
    P2avg_.resize(numbins_);
    for (int i=0;i<numbins_;i++)
    {
        Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 1.0/totalFrames);
        if (NumResPerBin_[i] < ignoreP2LessThan_)
        {
            Qtensor::matrix_mult_inplace(BinnedMatrix_[i], 0.0);
        }

        NumResPerBin_[i] /= totalFrames;
    }

    ASSERT((eigvec_.size() == BinnedMatrix_.size()), "failed.");
    for (int i=0;i<numbins_;i++)
    {
        auto ans = Qtensor::orderedeig_Qtensor(BinnedMatrix_[i]);
        P2avg_[i] = ans.second[0];

        for (int j=0;j<3;j++)
        {
            eigvec_[i][j] = ans.first[j][0];
        }
    }

    // finally, make the number of residues per bin to be 0
    for (int i=0;i<numbins_;i++)
    {
        if (NumResPerBin_[i] < ignoreP2LessThan_)
        {
            NumResPerBin_[i] = 0;
        }
    }
}