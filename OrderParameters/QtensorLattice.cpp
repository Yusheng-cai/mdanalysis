#include "QtensorLattice.hpp"

QtensorLattice::QtensorLattice(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadArrayNumber("LatticeShape", ParameterPack::KeyType::Required, lattice_shape_);
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    pack_.ReadNumber("cutoff", ParameterPack::KeyType::Required, cutoff_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailIndex_);
    headIndex_--;
    tailIndex_--;

    // initialize residue
    initializeResidueGroup(resname_);

    // resize the lattice 
    lattice_.resize(lattice_shape_);

    // cell grid
    cell_  = cellptr(new CellGrid(simstate_, cutoff_, 1));
}

void QtensorLattice::calculate()
{
    // assign the cell grid 
    const auto& res = getResidueGroup(resname_).getResidues();

    // resize uij and COM 
    uij_.resize(res.size());
    COM_.resize(res.size());


    // calculate the center of mass as well as the uij
    #pragma omp parallel for
    for (int i=0;i<COM_.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 distance;
        Real distsq;

        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, distsq);

        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    // calculate the cell indices 
    std::vector<std::vector<int>> cellIndices = cell_->calculateIndices(COM_);

    // fill the lattice
    #pragma omp parallel for 
    for (int i=0;i<lattice_.getSize();i++)
    {
        INT3 index3 = lattice_.getIndex3(i);
        Real3 pos;
        for (int j=0;j<3;j++)
        {
            pos[j] = index3[j] * dL_[j];
        }
        lattice_(index3) = pos;
    }

}

void QtensorLattice::update()
{
    // update the cell grid
    cell_->update();

    // update the sides 
    Real3 sides = simstate_.getSimulationBox().getSides();
    for (int i=0;i<3;i++)
    {
        dL_[i] = sides[i] / lattice_shape_[i];
    }
}