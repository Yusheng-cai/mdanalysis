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

        void CalculateDirector(std::vector<Real3>& uij, std::vector<Real>& norm, std::vector<Real>& indicators, Real& Ntilde, Real& N);
        void CalculateQtilde(const std::vector<Real3>& uij, Matrix& Qtensor, const std::vector<Real>& indicators, Real Ntilde);

        Real getN() const {return N_;}
        Real getNtilde() const {return Ntilde_;}

    protected:
        std::string tailgroupname_;
        std::size_t tailgroupsize_;

        std::string headgroupname_;
        std::string indicatorgroupname_;
        std::string pvname_;

        std::size_t headgroupsize_;

        Matrix Qtensor_;

        std::vector<Real3> uij_;
        std::vector<Real> norms_;
        std::vector<Real> indicators_;
        Real Ntilde_, N_;
};