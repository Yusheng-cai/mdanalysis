#include "OrderParameters.h"
#include "tools/Assert.h"
#include "Qtensor.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP_buffer.h"

#include <string>
#include <vector>

class P2cos: public OrderParameters
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        P2cos(const OrderParametersInput& input);
        virtual ~P2cos(){};
        virtual void calculate() override;
        virtual void update() override;
        std::pair<Real3,Real3> dP2cosdr(Real N, Real norm, Real3& eigvec, Real3& director);

        Real getP2cos() const {return P2cos_OP_;}

    private:
        std::string tailgroupname_;
        std::size_t tailgroupsize_;

        std::string headgroupname_;
        std::size_t headgroupsize_;

        Real P2cos_OP_;

        std::vector<Real3> uij_;
        std::vector<Real> norms_;

        std::array<Real,3> n_;
};