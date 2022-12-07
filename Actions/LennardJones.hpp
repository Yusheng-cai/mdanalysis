#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"

#include <cmath>

class LennardJones 
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        LennardJones(Real sigma, Real epsilon);

        Real calculate(const Real& r);

        Real3 calculate_force(const Real3& dist, const Real& r);

    private:
        Real sigma_;
        Real epsilon_;
        Real prefactor_;
};