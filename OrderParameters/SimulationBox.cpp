#include "SimulationBox.h"

SimulationBox::SimulationBox(Matrix& box)
{
    for (int i=0;i<3;i++)
    {
        length_[i] = box[i][i];
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
        center_[i] = length_[i] * 0.5;
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
        center_[i] = length_[i] * 0.5;
    }
}

SimulationBox::SimulationBox(Real3& length)
:length_(length)
{
    for (int i=0;i<3;i++)
    {
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
        center_[i] = length_[i] * 0.5;
    }
}

void SimulationBox::calculateDistance(const Real3& x1, const Real3& x2, Real3& distance, Real& sq_dist) const
{
    sq_dist = 0.0;
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

SimulationBox::Real3 SimulationBox::calculateShift(const Real3& x1, const Real3& ref) const
{
    Real3 shift;
    for (int i=0;i<3;i++)
    {
        Real dist = x1[i] - ref[i];

        if (dist < mhlength_[i]) { shift[i] = length_[i]; }
        else if (dist > hlength_[i]) { shift[i] = -length_[i];}
        else {shift[i] = 0.0;}
    }

    return shift;
}

SimulationBox::Real3 SimulationBox::calculateShift(const Real3& x1, const Real3& ref, const Real3& sides) const
{
    Real3 shift;
    for (int i=0;i<3;i++)
    {
        Real dist = x1[i] - ref[i];

        if (dist < -0.5 * sides[i]) { shift[i] = sides[i]; }
        else if (dist > 0.5 * sides[i]) { shift[i] = -sides[i];}
        else {shift[i] = 0.0;}
    }

    return shift;
}

SimulationBox::Real3 SimulationBox::shiftIntoBox(const Real3& x1) const 
{
    Real3 shift = calculateShift(x1, center_);

    Real3 ret = {};
    for (int i=0;i<3;i++)
    {
        ret[i] = x1[i] + shift[i];
    }

    return ret;
}

void SimulationBox::setBoxMatrix(const Matrix& box)
{
    // set the simulation matrix
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            box_[i][j] = box[i][j];
        }
    }

    // set the center of the box
    for (int i=0;i<3;i++)
    {
        center_[i] = box[i][i] * 0.5;
    }

    // set the length/half lengths of the box
    for (int i=0;i<3;i++)
    {
        length_[i] = box[i][i];
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}

void SimulationBox::setBoxMatrix(const Real3& sides)
{
    // set the simulation matrix
    box_ = {};
    for (int i=0;i<3;i++){
        box_[i][i] = sides[i];
    }

    // set the center of the box
    for (int i=0;i<3;i++){
        center_[i] = sides[i] * 0.5;
    }

    // set the length/half lengths of the box
    for (int i=0;i<3;i++){
        length_[i] = sides[i];
        hlength_[i] = 0.5*length_[i];
        mhlength_[i] = -hlength_[i];
    }
}

void SimulationBox::setCenter(const Real3& center){
    for (int i=0;i<3;i++){
        center_[i] = center[i];
    }
}