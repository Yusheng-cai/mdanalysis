#include "OrderParameters.h"
#include "tools/Assert.h"
#include "Qtensor.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP_buffer.h"
#include "liquid_crystal.h"
#include <cmath>

#include <string>
#include <vector>

class P2cos: public liquid_crystal
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        P2cos(const OrderParametersInput& input);
        virtual ~P2cos(){};

        virtual void calculate() override;
        virtual void update() override;
        std::pair<Real3,Real3> dP2cosdr(Real N, Real norm, Real3& director, Real3& n);

        Real getP2cos() const {return P2cos_OP_;}

    private:
        Real P2cos_OP_;

        std::array<Real,3> n_;
};