#include "QtensorCalc.h"

namespace CalculationRegistry
{
    registry_<QtensorCalc> pQtensorCalcRegister_("Qtensor");
}

QtensorCalc::QtensorCalc(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required, head_index_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required, tail_index_);
    head_index_--;
    tail_index_--;

    // read in the residue 
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_name_);
    initializeResidueGroup(residue_name_);

    // read in the array 
    pack_.ReadArrayNumber("array" , ParameterPack::KeyType::Optional, arr_);

    auto& res = getResidueGroup(residue_name_);
    size_  = res.size();
    uij_.resize(size_);

    // initialize Qtensor total --> Used to be averaged over the entire trajectory
    QtensorTot_.fill({});

    // register per iter output
    registerOutputFunction("Q", [this](std::string name) -> void {this -> printaverageQ(name);});
    registerPerIterOutputFunction("cos20", [this](std::ofstream& ofs) -> void {this -> printcos2thetaPerIter(ofs);});
    registerPerIterOutputFunction("cos0", [this](std::ofstream& ofs) -> void {this -> printcosPerIter(ofs);});
    registerPerIterOutputFunction("ev", [this](std::ofstream& ofs) -> void {this -> printevPerIter(ofs);});
    registerPerIterOutputFunction("p2", [this](std::ofstream& ofs) -> void {this -> printp2PerIter(ofs);});
    registerPerIterOutputFunction("qtensor", [this](std::ofstream& ofs) -> void {this -> printQtensorPerIter(ofs);});
    registerPerIterOutputFunction("cos2dir", [this](std::ofstream& ofs) -> void {this -> printcos2PerIter(ofs);});

    // register outputfile output
    registerOutputFileOutputs("biaxiality", [this](void) -> Real {return this -> getBiaxiality();});
    registerOutputFileOutputs("p2", [this](void) -> Real{return this -> getP2();});
    registerOutputFileOutputs("Qxx", [this](void) -> Real{return this -> getQxx();});
    registerOutputFileOutputs("Qxy", [this](void) -> Real {return this -> getQxy();});
    registerOutputFileOutputs("Qxz", [this](void) -> Real {return this -> getQxz();});
    registerOutputFileOutputs("Qyy", [this](void) -> Real {return this -> getQyy();});
    registerOutputFileOutputs("Qyz", [this](void) -> Real {return this -> getQyz();});
    registerOutputFileOutputs("v1x", [this](void) -> Real {return this -> getv1x();});
    registerOutputFileOutputs("v1y", [this](void) -> Real {return this -> getv1y();});
    registerOutputFileOutputs("v1z", [this](void) -> Real {return this -> getv1z();});

    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();
}

void QtensorCalc::printQtensorPerIter(std::ofstream& ofs)
{
    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);
    ofs << Qtensor_[0][0] << " " << Qtensor_[0][1] << " " << Qtensor_[0][2] << " " << Qtensor_[1][1] << " " << Qtensor_[1][2] \
    << " " << result.second[0] << " " << result.second[1] << " " << result.second[2] << "\n";
}

void QtensorCalc::printaverageQ(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    auto res = Qtensor::orderedeig_Qtensor(QtensorTot_);    
    Real3 eigenvector;
    for (int i=0;i<3;i++)
    {
        eigenvector[i] = res.first[i][0];
    }

    Real p2 = res.second[0];

    ofs << "# Qxx Qxy Qxz Qyy Qyz p2 v1x v1y v1z\n";
    ofs << QtensorTot_[0][0] << " " << QtensorTot_[0][1] << " " << QtensorTot_[0][2] << " " << QtensorTot_[1][1] << " " << QtensorTot_[1][2] << " " << p2 << \
    " " << eigenvector[0] << " " << eigenvector[1] << " " << eigenvector[2];

    ofs.close();
}

void QtensorCalc::printevPerIter(std::ofstream& ofs)
{
    Real3 v;
    for (int i=0;i<3;i++)
    {
        v[i] = eigenvector_[i][0];

        ofs << v[i] << " ";
    }

    ofs << "\n";
}

void QtensorCalc::printp2PerIter(std::ofstream& ofs)
{
    ofs << eigenval_[1] * (-2.0) << "\n";
}

void QtensorCalc::calculate()
{
    // zero out Qtensor 
    Qtensor_.fill({});

    auto& res = getResidueGroup(residue_name_).getResidues();
    int N = res.size();

    COM_.clear();
    COM_.resize(N);

    // First let's calculate the COM 
    #pragma omp parallel for 
    for (int i=0;i<N;i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
    }

    std::vector<int> InsideIndices = InsidePVIndices(COM_);

    int num = 0;
    for (int i=0;i<InsideIndices.size();i++)
    { 
        int index = InsideIndices[i];
        num ++;
        const auto& r = res[index];
        const auto& atoms = r.atoms_;
        Real3 distance;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(atoms[head_index_].positions_, atoms[tail_index_].positions_, distance, dist_sq);

        uij_[i] = Qtensor::normalize_director(distance);
        Matrix localQ = LinAlg3x3::LocalQtensor(uij_[i]);

        Qtensor::matrix_accum_inplace(Qtensor_, localQ);
    }

    // calculate the actualy Qtensor
    Qtensor::matrix_mult_inplace(Qtensor_, 1/(2.0*num));

    // Add the current Qtensor to Qtensortot
    Qtensor::matrix_accum_inplace(QtensorTot_, Qtensor_);

    // order from largest to smallest
    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);

    eigenvector_=  result.first;
    eigenval_   = result.second;

    p2_ = eigenval_[0];

    biaxiality_ = eigenval_[1] * 2.0 + eigenval_[0];
}

void QtensorCalc::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            QtensorTot_[i][j] /=  numframes;
        }
    }
}

void QtensorCalc::printcos2PerIter(std::ofstream& ofs)
{
    std::vector<Real> dottedval(uij_.size(),0.0);
    auto& res = getResidueGroup(residue_name_);

    for (int i=0;i<uij_.size();i++)
    {
        Real val = LinAlg3x3::DotProduct(arr_, uij_[i]);
        val = val * val;
        dottedval[i] = val;
    }

    int frameNum = simstate_.getFrameNumber();
    ofs << frameNum << " ";

    for (int i=0;i<dottedval.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            ofs << dottedval[i] << " ";
        }
    }

    ofs << "\n";
}

void QtensorCalc::printcosPerIter(std::ofstream& ofs)
{
    std::vector<Real> dottedval(uij_.size(),0.0);
    auto& res = getResidueGroup(residue_name_);

    for (int i=0;i<uij_.size();i++)
    {
        Real val = LinAlg3x3::DotProduct(arr_, uij_[i]);
        dottedval[i] = val;
    }

    int frameNum = simstate_.getFrameNumber();
    ofs << frameNum << " ";

    for (int i=0;i<dottedval.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            ofs << dottedval[i] << " ";
        }
    }

    ofs << "\n";
}

void QtensorCalc::printcos2thetaPerIter(std::ofstream& ofs)
{
    Real3 v;
    auto& res = getResidueGroup(residue_name_);
    for (int i=0;i<3;i++)
    {
        v[i] = eigenvector_[i][0];
    }

    std::vector<Real> dottedVal(uij_.size(),0.0);

    for (int i=0;i<uij_.size();i++)
    {
        Real val  = LinAlg3x3::DotProduct(v, uij_[i]);
        val = val*val;

        dottedVal[i] = val;
    }

    int frameNum = simstate_.getFrameNumber();
    ofs << frameNum << " ";

    for (int i=0;i<dottedVal.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            ofs << dottedVal[i] << " ";
        }
    }

    ofs << "\n";
}
