#pragma once

#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"

#include <vector>
#include <array>
#include <cmath>

class MorsePotential
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        MorsePotential(Real De, Real re, Real alpha);

        Real calculate(const Real& r);

        // return dU/dr
        Real3 calculate_force(const Real3& dist, const Real& r);

    private:
        Real De_;
        Real re_;
        Real alpha_;
};