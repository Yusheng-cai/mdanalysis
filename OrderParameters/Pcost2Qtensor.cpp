#include "Pcost2Qtensor.h"

namespace CalculationRegistry
{
    registry_<Pcost2Qtensor> pcost2QtensorRegister_("pcost2Qtensor");
}

Pcost2Qtensor::Pcost2Qtensor(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required, head_index_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required, tail_index_);
    head_index_--;
    tail_index_--;

    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_name_);
    addResidueGroup(residue_name_);

    auto& res = getResidueGroup(residue_name_);
    size_  = res.size();
    atomSize_ = res.getAtomSize();
    uij_.resize(size_);
    norm_.resize(size_,0.0);
    cost2Peratom_.resize(3, std::vector<Real>(atomSize_,0.0));
    cost2PerRes_.resize(3, std::vector<Real>(size_,0.0));

    auto bfPack = pack_.findParamPack("betafactor", ParameterPack::KeyType::Optional);
    if (bfPack != nullptr)
    {
        bf_ = bptr(new BetaFactorWriter(const_cast<ParameterPack&>(*bfPack)));
    }
    

    registerPerIterOutputFunction("cost20", [this](std::ofstream& ofs) -> void { this -> printcost20PerResIter(ofs);});
    registerPerIterOutputFunction("cost21", [this](std::ofstream& ofs) -> void { this -> printcost21PerResIter(ofs);});
    registerPerIterOutputFunction("cost22", [this](std::ofstream& ofs) -> void { this -> printcost22PerResIter(ofs);});
}

void Pcost2Qtensor::calculate()
{
    // zero out Qtensor 
    Qtensor_.fill({});

    auto& res = getResidueGroup(residue_name_);
    int N = res.size();

    for (int i=0;i<N;i++)
    {
        const auto& r = res[i];
        const auto& atoms = r.atoms_;
        Real3 distance;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(atoms[head_index_].positions_, atoms[tail_index_].positions_, distance, dist_sq);

        uij_[i] = Qtensor::normalize_director(distance);
        norm_[i] = std::sqrt(dist_sq);

        Matrix d = Qtensor::vec_dyadic(uij_[i], uij_[i]);
        d = Qtensor::matrix_sub(d, Qtensor::matrix_Identity());

        Qtensor::matrix_accum_inplace(Qtensor_, d);
    }

    Qtensor::matrix_mult_inplace(Qtensor_, 1/(2.0*N));

    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            eigenvector_[i][j] = result.first[j][i];
        }
    }

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<size_;j++)
        {
            Real cost2 = std::pow(Qtensor::vec_dot(uij_[j], eigenvector_[i]),2.0);
            cost2PerRes_[i][j] = cost2;
            int aSize = res[j].atoms_.size();

            for (int k=0;k<aSize;k++)
            {
                int atomIndex = res[j].atoms_[k].atomNumber_-1;
                cost2Peratom_[i][atomIndex] = cost2;
            }
        }
    }
}

void Pcost2Qtensor::printcost20PerResIter(std::ofstream& ofs)
{
    for (int i=0;i<cost2PerRes_[0].size();i++)
    {
        ofs << cost2PerRes_[0][i] << "\t";
    }

    ofs << "\n";
}

void Pcost2Qtensor::printcost21PerResIter(std::ofstream& ofs)
{
    for (int i=0;i<cost2PerRes_[1].size();i++)
    {
        ofs << cost2PerRes_[1][i] << "\t";
    }

    ofs << "\n";
}

void Pcost2Qtensor::printcost22PerResIter(std::ofstream& ofs)
{
    for (int i=0;i<cost2PerRes_[2].size();i++)
    {
        ofs << cost2PerRes_[2][i] << "\t";
    }

    ofs << "\n";
}

void Pcost2Qtensor::printOutputOnStep()
{
    if (bf_.get() != nullptr)
    {
        int timestep = simstate_.getFrameNumber();

        bf_ -> write(timestep, cost2Peratom_[0]);
    }
}