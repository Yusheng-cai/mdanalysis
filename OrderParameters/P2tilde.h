#pragma once
#include "OrderParameters.h"
#include "LinAlgTools.h"
#include "tools/CommonTypes.h"
#include "LiquidCrystal.h"

#include <vector>
#include <string>

class P2tilde: public LiquidCrystal
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        P2tilde(const OrderParametersInput& input);
        virtual ~P2tilde(){};

        virtual void calculate() override;

        Real getP2tilde() const {return p2tilde_;}
        Real getQxx() const {return Qtensor_[0][0];}
        Real getQxy() const {return Qtensor_[0][1];}
        Real getQxz() const {return Qtensor_[0][2];}
        Real getQyy() const {return Qtensor_[1][1];}
        Real getQyz() const {return Qtensor_[1][2];}
        Real geteig1() const {return eig1_;}
        Real geteig2() const {return eig2_;}
        Real getbiaxiality() const {return biaxiality_;}

    private:
        Real p2tilde_;
        Real3 v1_;
        Real biaxiality_;
        Real eig1_;
        Real eig2_;
};