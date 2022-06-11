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

void QtensorCalc::printaverageQ(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Qxx Qxy Qxz Qyy Qyz p2 v1x v1y v1z\n";
    ofs << QtensorTot_[0][0] << " " << QtensorTot_[0][1] << " " << QtensorTot_[0][2] << " " << QtensorTot_[1][1] << " " << QtensorTot_[1][2] << " " << p2_ << \
    " " << v0_[0] << " " << v0_[1] << " " << v0_[2];

    ofs.close();
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

        LinAlg3x3::normalize(distance);
        uij_[i] = distance;

        Matrix localQ = LinAlg3x3::LocalQtensor(uij_[i]);

        LinAlg3x3::matrix_accum_inplace(Qtensor_, localQ);
    }

    // calculate the actualy Qtensor
    LinAlg3x3::matrix_mult_inplace(Qtensor_, 1/(2.0*num));

    // Add the current Qtensor to Qtensortot
    LinAlg3x3::matrix_accum_inplace(QtensorTot_, Qtensor_);
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

    auto pair = LinAlg3x3::OrderEigenSolver(QtensorTot_);
    p2_ = pair.first[0];

    for (int i=0;i<3;i++)
    {
        v0_[i] = pair.second[i][0];

    }
}
