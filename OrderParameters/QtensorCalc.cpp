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

    // Read in the probeVolume
    pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvName_);

    // read in the residue 
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_name_);
    initializeResidueGroup(residue_name_);

    // read in the array 
    pack_.ReadArrayNumber("array" , ParameterPack::KeyType::Optional, arr_);

    auto& res = getResidueGroup(residue_name_);
    size_  = res.size();
    uij_.resize(size_);

    // register per iter output
    registerPerIterOutputFunction("cos20", [this](std::ofstream& ofs) -> void {this -> printcos2thetaPerIter(ofs);});
    registerPerIterOutputFunction("ev", [this](std::ofstream& ofs) -> void {this -> printevPerIter(ofs);});
    registerPerIterOutputFunction("p2", [this](std::ofstream& ofs) -> void {this -> printp2PerIter(ofs);});
    registerPerIterOutputFunction("qtensor", [this](std::ofstream& ofs) -> void {this -> printQtensorPerIter(ofs);});
    registerPerIterOutputFunction("cos2dir", [this](std::ofstream& ofs) -> void {this -> printcos2PerIter(ofs);});
}

void QtensorCalc::printQtensorPerIter(std::ofstream& ofs)
{
    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);
    ofs << Qtensor_[0][0] << " " << Qtensor_[0][1] << " " << Qtensor_[0][2] << " " << Qtensor_[1][1] << " " << Qtensor_[1][2] \
    << " " << result.second[0] << " " << result.second[1] << " " << result.second[2] << "\n";
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

    auto& pv  = simstate_.getProbeVolume(pvName_);
    COM_.clear();
    COM_.resize(N);

    // First let's calculate the COM 
    #pragma omp parallel for 
    for (int i=0;i<N;i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
    }

    int num = 0;
    for (int i=0;i<N;i++)
    { 
        auto result = pv.calculate(COM_[i]);

        if (result.hx_ == 1)
        {
            num ++;
            const auto& r = res[i];
            const auto& atoms = r.atoms_;
            Real3 distance;
            Real dist_sq;

            simstate_.getSimulationBox().calculateDistance(atoms[head_index_].positions_, atoms[tail_index_].positions_, distance, dist_sq);

            uij_[i] = Qtensor::normalize_director(distance);
            Matrix localQ = LinAlg3x3::LocalQtensor(uij_[i]);

            Qtensor::matrix_accum_inplace(Qtensor_, localQ);
        }
    }

    // calculate the actualy Qtensor
    Qtensor::matrix_mult_inplace(Qtensor_, 1/(2.0*num));

    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);

    eigenvector_=  result.first;
    eigenval_   = result.second;
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
