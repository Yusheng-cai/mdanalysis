#include "SimulationBox.h"

SimulationBox::SimulationBox(Matrix& box)
{
    for (int i=0;i<3;i++)
    {
        length_[i] = box[i][i];
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}

SimulationBox::SimulationBox(Real lx, Real ly, Real lz)
{
    length_[0] = lx;
    length_[1] = ly;
    length_[2] = lz;

    for (int i=0;i<3;i++)
    {
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}

SimulationBox::SimulationBox(Real3& length)
:length_(length)
{
    for (int i=0;i<3;i++)
    {
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}

void SimulationBox::calculateDistance(const Real3& x1, const Real3& x2, Real3& distance, Real& sq_dist)
{
    for (int i=0;i<3;i++)
    {
        Real dist = x1[i] - x2[i];
        if (dist > hlength_[i]) {dist -= length_[i];}
        else if (dist < mhlength_[i]) {dist += length_[i];}
        else { dist = dist;}

        distance[i] = dist;
        sq_dist += dist*dist;
    }
}

SimulationBox::Real3 SimulationBox::calculateShift(const Real3& x1, const Real3& ref)
{
    Real3 shift;
    for (int i=0;i<3;i++)
    {
        Real dist = x1[i] - ref[i];

        if (dist < mhlength_[i]) { shift[i] = length_[i]; }
        else if (dist > mhlength_[i]) { shift[i] = -length_[i];}
        else {shift[i] = 0.0;}
    }

    return shift;
}

void SimulationBox::setBoxMatirx(const Matrix& box)
{
    for (int i=0;i<3;i++)
    {
        length_[i] = box[i][i];
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}