#include "OrderParameters.h"
#include "tools/Assert.h"
#include "Qtensor.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP_buffer.h"

#include <string>
#include <vector>

class liquid_crystal: public OrderParameters
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        liquid_crystal(const OrderParametersInput& input);
        virtual ~liquid_crystal(){};

        void getUij();

    private:
        std::string tailgroupname_;
        std::size_t tailgroupsize_;

        std::string headgroupname_;
        std::size_t headgroupsize_;

        Real P2_OP_;
        Real3 v1_;

        std::vector<Real3> uij_;
        std::vector<Real> norms_;

        OpenMP::OpenMP_buffer<Matrix> Qtensor_buffer_;
};