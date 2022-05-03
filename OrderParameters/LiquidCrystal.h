#pragma once

#include "OrderParameters.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "parallel/OpenMP_buffer.h"

#include <string>
#include <vector>

class LiquidCrystal: public OrderParameters
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        LiquidCrystal(const OrderParametersInput& input);
        virtual ~LiquidCrystal(){};

        void getUij();
        void calcQtensor();

    protected:
        std::string tailgroupname_;
        std::size_t tailgroupsize_;

        std::string headgroupname_;
        std::size_t headgroupsize_;

        Matrix Qtensor_;

        std::vector<Real3> uij_;
        std::vector<Real> norms_;
};