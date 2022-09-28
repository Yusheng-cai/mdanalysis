#include "SmecticOrder.h"

SmecticOrder::SmecticOrder(const CalculationInput& input)
: Calculation(input)
{
    // initialize the residue group
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residuegroup_);
    initializeResidueGroup(residuegroup_);

    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;
}

void SmecticOrder::calculate()
{
    const auto& resgroup = getResidueGroup(residuegroup_).getResidues();
    uij_.clear();
    uij_.resize(resgroup.size());

    #pragma omp parallel for
    for (int i=0;i<resgroup.size();i++)
    {
        COM_[i] = calcCOM(resgroup[i]);
        Real3 headpos = resgroup[i].atoms_[headindex_].positions_;
        Real3 tailpos = resgroup[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distsq);

        distance = distance / distsq;
        uij_[i] = distance;
    }

    // calculate the Qtensor
    Matrix Q;
    Q.fill({});
    for (int i=0;i<resgroup.size();i++)
    {
        Matrix Qlocal = LinAlg3x3::LocalQtensor(uij_[i]);
        Q = Q + Qlocal;
    }
}

void SmecticOrder::update()
{

}