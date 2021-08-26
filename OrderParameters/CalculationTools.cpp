#include "CalculationTools.h"

std::vector<CalculationTools::Real3> CalculationTools::getCOM(const std::vector<Molecule::residue>& residueGroup, const SimulationState& simstate)
{
    std::vector<Real3> COM_;
    COM_.resize(residueGroup.size());

    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    #pragma omp parallel for
    for (int i=0;i<residueGroup.size();i++)
    {
        auto& pos1 = residueGroup[i].atoms_[0].positions_;
        Real massTot = 0.0;
        Real3 COM_pos = {{0,0,0}};
        for (int j=1;j<residueGroup[i].atoms_.size();j++)
        {
            Real3 distance;
            Real distsq;
            simbox.calculateDistance(residueGroup[i].atoms_[j].positions_, pos1, distance, distsq);

            for (int k=0;k<3;k++)
            {
                COM_pos[k] += distance[k] * residueGroup[i].atoms_[j].mass_;
            }

            massTot += residueGroup[i].atoms_[j].mass_;
        }

        for (int j=0;j<3;j++)
        {
            COM_pos[j] = COM_pos[j]/massTot + pos1[j];
        }

        // check if it is outside of the box, if it is, then move it inside the box
        Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());

        for (int j=0;j<3;j++)
        {
            COM_pos[j] += shift[j];
        }

        COM_[i] = COM_pos;
    }

    return COM_;
}