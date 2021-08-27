#include "CalculationTools.h"

CalculationTools::Real3 CalculationTools::getCOM(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices_)
{
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices_[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real massTot = 0;
    Real3 COM_pos = {{0,0,0}};
    

    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++)
    {
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos;

        for (int k=0;k<3;k++)
        {
            shiftedPos[k] = residueGroup.atoms_[index].positions_[k] + shiftWRTatom1[k];
        }

        for (int k=0;k<3;k++)
        {
            COM_pos[k] += shiftedPos[k] * residueGroup.atoms_[index].mass_;
        }

        massTot += residueGroup.atoms_[index].mass_;
    }

    for (int j=0;j<3;j++)
    {
        COM_pos[j] = COM_pos[j]/massTot;
    }

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());

    for (int j=0;j<3;j++)
    {
        COM_pos[j] += shift[j];
    }

    return COM_pos;
}