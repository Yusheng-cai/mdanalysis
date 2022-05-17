#include "P2tilde.h"

namespace OrderParametersRegistry
{
    registry_<P2tilde> registerP2tilde("p2tilde");
}

P2tilde::P2tilde(const OrderParametersInput& input)
:LiquidCrystal(input)
{
    registerOutput("p2tilde", [this](void)->Real {return this->getP2tilde();});

    registerOutput("qxx", [this](void) -> Real {return this->getQxx();});
    registerOutput("qxy", [this](void) -> Real {return this->getQxy();});
    registerOutput("qxz", [this](void) -> Real {return this -> getQxz();});
    registerOutput("qyy", [this](void)-> Real {return this -> getQyy();});
    registerOutput("qyz", [this](void) -> Real {return this -> getQyz();});
    registerOutput("biaxiality", [this](void)->Real {return this->getbiaxiality();});
    registerOutput("eig1", [this](void) -> Real {return this->geteig1();});
    registerOutput("eig2", [this](void) -> Real {return this->geteig2();});
}

void P2tilde::calculate()
{
    // calculate director 
    CalculateDirector(uij_, norms_, indicators_, Ntilde_, N_);

    // calculate Qtensor 
    CalculateQtilde(uij_, Qtensor_, indicators_, Ntilde_);

    // solve the Qtensor 
    auto result = LinAlg3x3::OrderEigenSolver(Qtensor_);

    // save the results
    p2tilde_ = result.first[0]; 
    eig1_    = result.first[1];
    eig2_    = result.first[2];
    biaxiality_ = eig1_ * 2.0 + p2tilde_;

    v1_ = result.second[0];
}