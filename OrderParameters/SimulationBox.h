#pragma once
#include "tools/CommonTypes.h"

// currently only support rectangular box
class SimulationBox
{
    public:
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;
        using Matrix = CommonTypes::Matrix;

        SimulationBox();
        SimulationBox(Matrix& box);
        SimulationBox(Real3& length);
        SimulationBox(Real lx, Real ly, Real lz);

        // setter 
        void setBoxMatirx(const Matrix& box);

        // calculate pbc corrected distance between x1 and x2
        void calculateDistance(const Real3& x1, const Real3& x2, Real3& distance, Real& sq_dist);

        // calculate the shift between x1 and x2
        Real3 calculateShift(const Real3& x1, const Real3& ref);

        // setters
        void setSimulationBox(const Matrix& box){box_ = box;}

    private:
        Matrix box_;

        // Length of the box (lx, ly, lz)
        Real3 length_;

        // half of the length of the box (lx/2, ly/2, lz/2)
        Real3 hlength_;

        // minus half of the length of the box 
        Real3 mhlength_;
};