#include "OrderParameters.h"
#include "tools/Assert.h"
#include "Qtensor.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP_buffer.h"
#include "liquid_crystal.h"

#include <string>
#include <vector>

class P2: public liquid_crystal
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        P2(const OrderParametersInput& input);
        virtual ~P2(){};
        virtual void calculate() override;
        virtual void update() override;
        std::pair<Real3,Real3> dP2dr(Real N, Real norm, Real3& eigvec, Real3& director);

        Real getP2() const {return P2_OP_;}
        Real getv1x() const {return v1_[0];}
        Real getv1y() const {return v1_[1];}
        Real getv1z() const {return v1_[2];}
        Real getQxx() const {return Qtensor_[0][0];}
        Real getQxy() const {return Qtensor_[0][1];}
        Real getQxz() const {return Qtensor_[0][2];}
        Real getQyy() const {return Qtensor_[1][1];}
        Real getQyz() const {return Qtensor_[1][2];}

    private:
        Real P2_OP_;
        Real3 v1_;
};