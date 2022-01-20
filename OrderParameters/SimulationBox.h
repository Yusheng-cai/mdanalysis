#pragma once
#include "tools/CommonTypes.h"

// currently only support rectangular box
class SimulationBox
{
    public:
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;
        using Matrix = CommonTypes::Matrix;

        SimulationBox(){};
        SimulationBox(Matrix& box);
        SimulationBox(Real3& length);
        SimulationBox(Real lx, Real ly, Real lz);
        ~SimulationBox(){};

        // setter 
        void setBoxMatrix(const Matrix& box);

        // calculate pbc corrected distance between x1 and x2
        void calculateDistance(const Real3& x1, const Real3& x2, Real3& distance, Real& sq_dist) const;

        // calculate the shift between x1 and x2
        Real3 calculateShift(const Real3& x1, const Real3& ref) const;

        // shift the position into box
        Real3 shiftIntoBox(const Real3& x1) const;

        // get simulation box
        const Matrix& getBox() const {return box_;}

        // get the volume of the box
        Real getVolume() const
        {
            Real3 sides = getSides();

            return sides[0]*sides[1]*sides[2];
        }

        // get the 3 sides of the cubic box
        Real3 getSides() const {return length_;}

        // get the center of the box
        Real3 getCenter() const {return center_;}

    private:
        Matrix box_;

        // Length of the box (lx, ly, lz)
        Real3 length_;

        // half of the length of the box (lx/2, ly/2, lz/2)
        Real3 hlength_;

        // minus half of the length of the box 
        Real3 mhlength_;

        // center of the box
        Real3 center_;
};