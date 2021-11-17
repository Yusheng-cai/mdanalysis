#include "LocalOrder.h"

namespace CalculationRegistry
{
    registry_<LocalOrder> registerLocalOrder("localorder");
}

LocalOrder::LocalOrder(const CalculationInput& input)
: Calculation(input)
{
    registerPerIterOutputFunction("p2", [this](std::ofstream& ofs) -> void {this -> printLocalOrderBetaFactor(ofs);});
    registerPerIterOutputFunction("biaxiality", [this](std::ofstream& ofs) -> void {this -> printLocalBiaxialityBetaFactor(ofs);});
    registerPerIterOutputFunction("director", [this](std::ofstream& ofs) -> void {this -> printLocaldirector(ofs);});

    // read in the radius 
    input.pack_.ReadNumber("radius", ParameterPack::KeyType::Required, radius_);

    // read in the residue group
    input.pack_.ReadString("residue", ParameterPack::KeyType::Optional, residueName_);
    initializeResidueGroup(residueName_);

    // read in the head & tail index of the residue 
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required, headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required, tailIndex_);
    input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    headIndex_--;
    tailIndex_--;
}

void LocalOrder::calculate()
{
    // get all the residues
    auto& res = getResidueGroup(residueName_).getResidues();
    COM_.resize(res.size());
    uij_.resize(res.size());
    localP2_.resize(res.size());
    localBiaxiality_.resize(res.size());
    localdirector_.resize(res.size());

    // obtain the COM 
    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, distsq);

        LinAlg3x3::normalize(distance);
        uij_[i] = distance;
    }

    NeighborIndicesBuffer_.set_master_object(neighborIndices_);
    #pragma omp parallel
    {
        auto& localbuffer = NeighborIndicesBuffer_.access_buffer_by_id();
        localbuffer.clear();
        localbuffer.resize(res.size());

        #pragma omp for 
        for (int i=0;i<COM_.size();i++)
        {
            for (int j=i+1;j<COM_.size();j++)
            {
                // find the distance between ith and jth COM 
                Real3 distance;
                Real val;
                simstate_.getSimulationBox().calculateDistance(COM_[i], COM_[j], distance, val);

                if (std::sqrt(val) <= radius_)
                {
                    localbuffer[i].push_back(j);
                    localbuffer[j].push_back(i);
                }
            }
        }
    }

    for (int i=0;i<res.size();i++)
    {
        for ( auto it = NeighborIndicesBuffer_.beginworker(); it != NeighborIndicesBuffer_.endworker(); it ++)
        {
            neighborIndices_[i].insert(neighborIndices_[i].end(), (*it)[i].begin(), (*it)[i].end());
        }
    }

    ASSERT((neighborIndices_.size() == res.size()), "There is a size mismatch.");

    // calculate the local Qtensors 
    #pragma omp parallel for 
    for (int i=0;i<COM_.size();i++)
    {
        Matrix QtensorLocal;
        QtensorLocal.fill({});

        Matrix Q = LinAlg3x3::LocalQtensor(uij_[i]);
        Qtensor::matrix_accum_inplace(QtensorLocal, Q);

        // iterate over the neighbors 
        for (int j=0;j<neighborIndices_[i].size();j++)
        {
            int index = neighborIndices_[i][j];

            Matrix Q  = LinAlg3x3::LocalQtensor(uij_[index]);

            Qtensor::matrix_accum_inplace(QtensorLocal, Q);
        }

        Qtensor::matrix_mult_inplace(QtensorLocal, 0.5/(neighborIndices_[i].size() + 1));

        // find the eigenvalue and eigenvector 
        auto res = Qtensor::orderedeig_Qtensor(QtensorLocal);
        localP2_[i] = res.second[0];
        localBiaxiality_[i] = res.second[1] * 2.0 + res.second[0]; 

        // get the local director dotted with the corresponding user provided direction
        Real dotProduct = LinAlg3x3::DotProduct(res.first[0], arr_);
        localdirector_[i] = dotProduct * dotProduct;
    }
}

void LocalOrder::printLocaldirector(std::ofstream& ofs)
{
    auto& res = getResidueGroup(residueName_);
    std::vector<Real> localdirectorPeratom_(res.getAtomSize(),0.0);

    int index = 0;
    for (int i=0;i<res.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            localdirectorPeratom_[index] = localdirector_[i];
            index ++;
        }
    }

    int timeframe = simstate_.getFrameNumber();

    ofs << timeframe << " ";

    for (int i=0;i<localdirectorPeratom_.size();i++)
    {
        ofs << localdirectorPeratom_[i] << " ";
    }

    ofs << "\n";

}

void LocalOrder::printLocalOrderBetaFactor(std::ofstream& ofs)
{
    auto& res = getResidueGroup(residueName_);
    std::vector<Real> localP2Peratom_(res.getAtomSize(),0.0);

    int index = 0;
    for (int i=0;i<res.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            localP2Peratom_[index] = localP2_[i];
            index ++;
        }
    }

    int timeframe = simstate_.getFrameNumber();

    ofs << timeframe << " ";

    for (int i=0;i<localP2Peratom_.size();i++)
    {
        ofs << localP2Peratom_[i] << " ";
    }

    ofs << "\n";
}

void LocalOrder::printLocalBiaxialityBetaFactor(std::ofstream& ofs)
{
    auto& res = getResidueGroup(residueName_);
    std::vector<Real> localBiaxialityPeratom_(res.getAtomSize(),0.0);

    int index = 0;
    for (int i=0;i<res.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            localBiaxialityPeratom_[index] = localBiaxiality_[i]; 
            index ++;
        }
    }

    int timeframe = simstate_.getFrameNumber();

    ofs << timeframe << " ";

    for (int i=0;i<localBiaxialityPeratom_.size();i++)
    {
        ofs << localBiaxialityPeratom_[i] << " ";
    }

    ofs << "\n";
}