#include "Frame.h"

void Frame::setNumAtoms(int num_atoms)
{
    ASSERT((num_atoms != 0), "The number of the atoms equals to 0 and is not valid.");

    num_atoms_ = num_atoms;
    positions_.resize(num_atoms);
    velocities_.resize(num_atoms);
    forces_.resize(num_atoms);
}