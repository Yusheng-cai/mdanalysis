#include "SlabQtensor.h"

namespace CalculationRegistry
{
    registry_<SlabQtensor> registerSlabQtensor("SlabQtensor");
}

SlabQtensor::SlabQtensor(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailIndex_);

    registerOutputFunction("p2z", [this](std::string name) -> void{this -> printP2z(name);});
    registerOutputFunction("evbeta", [this](std::string name) -> void {this -> printevBeta(name);});
    registerPerIterOutputFunction("QtensorZ", [this](std::ofstream& ofs) -> void {this -> printQtensorZPerIter(ofs);});

    // add the residue group to the system
    initializeResidueGroup(residueName_);
    headIndex_--;
    tailIndex_--;
    
    // initialize bin
    auto binPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);
    // create the bin object , if specified, then used user specified bin, else we use min max of the COM
    if (binPack != nullptr){
        bin_ = Binptr(new Bin(*binPack));
        numbins_ = bin_ -> getNumbins();
    }
    else{
        input.pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numbins_);
        input.pack_.ReadNumber("above", ParameterPack::KeyType::Optional, above_);
        input.pack_.ReadNumber("below", ParameterPack::KeyType::Optional, below_);
        bin_ = Binptr(new Bin());
        usingMinMaxBins_ = true;
        BinLocation_.resize(numbins_, 0.0);
    }

    // initialize direction
    auto it = MapNameToDirection.find(direction_);
    ASSERT((it != MapNameToDirection.end()), "The direction " << direction_ << " is not recognized.");
    index_ = it -> second;
    
    // resize the binned matrix to number of bins
    // make all the matrix in vector zero
    QtensorZ_.resize(numbins_, {});
    P2_perIter_.resize(numbins_, 0);
    eigvec_perIter_.resize(numbins_, {});
    ResiduePerBin_perIter_.resize(numbins_,0.0);

    // resize number of residues per bin
    NumResPerBin_.resize(numbins_,0.0);

    // res index to bin index 
    auto& res = getResidueGroup(residueName_);
    ResIndexToBinIndex_.resize(res.size(),0);
    BetaFactors_.resize(res.getAtomSize(),0.0);

    // shift the box 
    input.pack_.ReadArrayNumber("box_shift", ParameterPack::KeyType::Optional, box_shift_);
}

void SlabQtensor::printP2zbeta(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    int numframes = simstate_.getTotalFrames();
    auto& res = getResidueGroup(residueName_).getResidues();
    int totalatoms = getResidueGroup(residueName_).getAtomSize();

    std::vector<Real> betaFactors(totalatoms,0.0);

    for (int j=0;j<res.size();j++){
        int binindex = ResIndexToBinIndex_[j];

        Real val = P2_[binindex];

        for (int k=0;k<res[j].atoms_.size();k++){
            int index = res[j].atoms_[k].atomNumber_ - 1;
            betaFactors[index] = val;
        }
    }
}

void SlabQtensor::printevBeta(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    auto& res = getResidueGroup(residueName_).getResidues();
    int totalatoms = getResidueGroup(residueName_).getAtomSize();

    std::vector<Real> betaFactors(totalatoms,0.0);

    for (int j=0;j<res.size();j++){
        int binindex = ResIndexToBinIndex_[j];

        Real val = eigvec_[binindex][2] * eigvec_[binindex][2];

        for (int k=0;k<res[j].atoms_.size();k++){
            int index = res[j].atoms_[k].atomNumber_ - 1;
            betaFactors[index] = val;
        }
    }

    ofs << 0  << " ";

    for (int j=0;j<betaFactors.size();j++){
        ofs << betaFactors[j] << " ";
    }
    ofs << "\n";

    ofs.close();
}

void SlabQtensor::binUsingMinMax()
{
    Real slight_shift=1e-3;

    std::vector<Real> zdir;
    for (int i=0;i<COM_.size();i++){
        if ((COM_[i][index_] > above_) && (COM_[i][index_] < below_)){
            zdir.push_back(COM_[i][index_]);
        }
    }

    Real min   = Algorithm::min(zdir) - slight_shift;
    Real max   = Algorithm::max(zdir) + slight_shift;

    Range range = {{min, max}};

    bin_->update(range, numbins_);

    // sum over the bin location 
    for (int i=0;i<numbins_;i++){
        BinLocation_[i] += bin_ -> getCenterLocationOfBin(i);
    }
    std::cout << "Min = " << min << ", Max = " << max << "\n";
}

void SlabQtensor::update(){
    Matrix curr_box    = simstate_.getSimulationBox().getBox();
    Real3 curr_center = simstate_.getSimulationBox().getCenter();
    Real3 update_center = curr_center + box_shift_;
    shifted_box_.setBoxMatrix(curr_box);
    shifted_box_.setCenter(update_center);
}

void SlabQtensor::calculate()
{
    // get Qtensor Z for this iteration
    std::vector<Matrix> QtensorZ_Iter;
    QtensorZ_Iter.resize(numbins_, {});
    std::fill(ResiduePerBin_perIter_.begin(), ResiduePerBin_perIter_.end(), 0.0);

    // obtain the residue group by its name
    const auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());

    // obtain the center of mass of each of the residues
    #pragma omp parallel for 
    for (int i=0;i<res.size();i++){
        COM_[i] = shifted_box_.shiftIntoBox(CalculationTools::getCOM(res[i], simstate_, COMIndices_));
    }

    // update the bins using COM min and max 
    if (usingMinMaxBins_){
        binUsingMinMax();
    }

    // Bin the COMs
    for (int i=0;i<COM_.size();i++){
        // get where the COM is --> index (x,y,z)
        Real num = COM_[i][index_];

        // only perform these operations of num is in range of bin
        if (bin_->isInRange(num)){
            int binNum = bin_ -> findBin(num);

            // calculate the uij 
            Real3 diff;
            Real diff_sq;
            Real3 headPos = res[i].atoms_[headIndex_].positions_;
            Real3 tailPos = res[i].atoms_[tailIndex_].positions_;
            simstate_.getSimulationBox().calculateDistance(headPos, tailPos, diff, diff_sq);

            ResIndexToBinIndex_[i] = binNum;
            
            // normalize the director and calculate the dyad
            LinAlg3x3::normalize(diff);
            Matrix Qlocal  = LinAlg3x3::LocalQtensor(diff);

            // accumulate matrix 
            LinAlg3x3::matrix_accum_inplace(QtensorZ_Iter[binNum], Qlocal);

            // accumulate the number of residue per bin 
            ResiduePerBin_perIter_[binNum] += 1;
            NumResPerBin_[binNum] += 1;
        }
    }

    // sum up the matrices from per iter to the total one
    for (int i=0;i<QtensorZ_Iter.size();i++){
        LinAlg3x3::matrix_accum_inplace(QtensorZ_[i], QtensorZ_Iter[i]);
    }

    // get the eigenvalue and eigenvector
    for (int i=0;i<QtensorZ_Iter.size();i++){
        // 0.5 * Qtensor 
        if (ResiduePerBin_perIter_[i] > 0){
            LinAlg3x3::matrix_mult_inplace(QtensorZ_Iter[i], 1.0/(2.0 * ResiduePerBin_perIter_[i]));

            // find the eigenvalue and eigenvector
            auto res = LinAlg3x3::OrderEigenSolver(QtensorZ_Iter[i]);

            // put the eigenvalue and eigenvector to vector
            P2_perIter_[i] = res.first[0];
            Real3 eigv;
            for (int j=0;j<3;j++){eigv[j] = res.second[j][0];}
            eigvec_perIter_[i] = eigv;
        }
    }
}

void SlabQtensor::printP2z(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    ofs << "#z\tp2\tNumber\tQxx\tQxy\tQxz\tQyy\tQyz\tQzz\tnx\tny\tnz\n";

    for (int i=0;i<numbins_;i++){
        Real binlocation;
        if (usingMinMaxBins_){
            binlocation = BinLocation_[i];
        }
        else{
            binlocation = bin_->getCenterLocationOfBin(i);
        }

        ofs << binlocation << "\t" << P2avg_[i] << "\t" << NumResPerBin_[i];
        ofs << "\t" << QtensorZ_[i][0][0] << "\t" << QtensorZ_[i][0][1] << "\t" << QtensorZ_[i][0][2]; 
        ofs << "\t" << QtensorZ_[i][1][1] << "\t" << QtensorZ_[i][1][2];
        ofs << "\t" << QtensorZ_[i][2][2];
        ofs << "\t" << eigvec_[i][0] << "\t" << eigvec_[i][1] << "\t" << eigvec_[i][2] << "\n";
    }
    ofs.close();
}

void SlabQtensor::printQtensorZPerIter(std::ofstream& ofs){
    for (int i=0;i<ResiduePerBin_perIter_.size();i++){
        ofs << ResiduePerBin_perIter_[i] << " ";
    }

    for (int i=0;i<P2_perIter_.size();i++){
        ofs << P2_perIter_[i] << " ";
    }

    for (int i=0;i<eigvec_perIter_.size();i++){
        for (int j=0;j<3;j++){
            ofs << eigvec_perIter_[i][j] << " ";
        }
    }

    ofs << "\n";
}

void SlabQtensor::printp2zPerIter(std::ofstream& ofs){
    for (int i=0;i<ResiduePerBin_perIter_.size();i++){
        ofs << ResiduePerBin_perIter_[i] << " ";
    }

    for (int i=0;i<P2_perIter_.size();i++){
        ofs << P2_perIter_[i] << " ";
    }

    ofs << "\n";
}

void SlabQtensor::finishCalculate()
{
    int totalFrames = simstate_.getTotalFrames();

    P2avg_.resize(numbins_);
    eigvec_.resize(numbins_);

    // average the Qtensor over time
    for (int i=0;i<numbins_;i++){
        if (NumResPerBin_[i] != 0){
            LinAlg3x3::matrix_mult_inplace(QtensorZ_[i], 1.0/(2.0*NumResPerBin_[i]));
        }

        // obtain normalized number residues per bin
        NumResPerBin_[i] /= totalFrames;

        // average over the bin location
        if (usingMinMaxBins_){
            BinLocation_[i] = BinLocation_[i] / totalFrames;
        }
    }

    // find the eigenvector and eigenvalue of the Qtensor
    for (int i=0;i<numbins_;i++){
        if (NumResPerBin_[i] != 0){
            auto ans = LinAlg3x3::OrderEigenSolver(QtensorZ_[i]);
            P2avg_[i] = ans.first[0];

            for (int j=0;j<3;j++){
                eigvec_[i][j] = ans.second[j][0];
            }
        }
    }
}