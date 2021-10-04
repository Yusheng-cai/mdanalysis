#include "P2.h"

namespace OrderParametersRegistry
{
    registry_<P2> registerP2("p2");
}

P2::P2(const OrderParametersInput& input)
:liquid_crystal(input)
{
    registerOutput("p2",[this](void)-> Real {return this->getP2();});
    registerOutput("p20", [this](void) -> Real {return this->getP20();});
    registerOutput("p21", [this](void) -> Real {return this->getP21();});
    registerOutput("p22", [this](void) -> Real {return this->getP22();});
    registerOutput("qxx", [this](void)->Real {return this->getQxx();});
    registerOutput("qxy", [this](void)->Real {return this->getQxy();});
    registerOutput("qxz", [this](void)->Real {return this ->getQxz();});
    registerOutput("qyy", [this](void)->Real {return this->getQyy();});
    registerOutput("qyz", [this](void) ->Real {return this->getQyz();});
    registerOutput("v0x", [this](void) -> Real {return this -> getv0x();});
    registerOutput("v0y", [this](void) -> Real {return this -> getv0y();});
    registerOutput("v0z", [this](void) -> Real {return this -> getv0z();});
}

void P2::calculate()
{
    getUij();

    calcQtensor();

    // find the p2 variable as well as the eigenvalue
    auto orderedEigPair = Qtensor::orderedeig_Qtensor(Qtensor_);

    P2_OP_ = -2.0*orderedEigPair.second[1];
    for (int i=0;i<3;i++)
    {
        v0_[i]    = orderedEigPair.first[i][0];
        v1_[i]    = orderedEigPair.first[i][1];
        v2_[i]    = orderedEigPair.first[i][2];
    }

    // calculate p2 in each of the directions of the eigenvectors 
    p2_0_ = 0.0;
    p2_1_ = 0.0;
    p2_2_ = 0.0;
    for (int i=0;i<uij_.size();i++)
    {
        p2_0_ += 3*std::pow(Qtensor::vec_dot(uij_[i], v0_),2.0) - 1;
        p2_1_ += 3*std::pow(Qtensor::vec_dot(uij_[i], v1_),2.0) - 1;
        p2_2_ += 3*std::pow(Qtensor::vec_dot(uij_[i], v2_),2.0) - 1;
    }

    auto& headatomgroup_ = getAtomGroup(headgroupname_);
    auto& tailatomgroup_ = getAtomGroup(tailgroupname_);

    auto& headatoms_ = headatomgroup_.getAtoms();
    auto& tailatoms_ = tailatomgroup_.getAtoms();

    Real N = tailgroupsize_;

    p2_0_ /= (2.0*N); 
    p2_1_ /= (2.0*N);
    p2_2_ /= (2.0*N);

    auto headatomDerivatives_ = accessDerivatives(headgroupname_);
    auto tailatomDerivatives_ = accessDerivatives(tailgroupname_);

    // clear all the derivatives stored in derivativesOutputs
    clearDerivativesOutputs();

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<tailgroupsize_;i++)
        {
            Real3 uij = uij_[i];
            Real norm = norms_[i];
            auto derivative = dP2dr(N, norm, v1_, uij);

            auto headatom = headatoms_[i];
            auto tailatom = tailatoms_[i];

            headatomDerivatives_.insertOMP(headatom.index, derivative.first);
            tailatomDerivatives_.insertOMP(tailatom.index, derivative.second);
        } 
    }
    headatomDerivatives_.CombineAndClearOMPBuffer();
    tailatomDerivatives_.CombineAndClearOMPBuffer();

    // for (int i=0;i<headatomDerivatives_.size();i++)
    // {
    //     auto deriv = headatomDerivatives_.getAtomDerivativeByIndex(i);

    //     std::cout << "For atom index " << deriv.index << ":";
    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << deriv.derivatives[j] << " ";
    //     }
    //     std::cout << "\n";
    // }
}

std::pair<P2::Real3,P2::Real3> P2::dP2dr(Real N, Real norm, Real3& eigvec, Real3& director)
{
    std::pair<Real3,Real3> pair_;

    // (v1 \dot u)
    Real dotproduct = Qtensor::vec_dot(eigvec, director);
    Real factor = - 6.0/(N*norm);

    Real3 output;
    output.fill(0);

    for (int i=0;i<3;i++)
    {
        output[i] += eigvec[i] - director[i]*dotproduct;
    }

    output = Qtensor::vec_mult(factor, output);
    Real3 output2 = Qtensor::vec_mult(-1.0, output);

    pair_.first = output;
    pair_.second = output2;

    return pair_;
}

void P2::update()
{
    P2_OP_ = 0;
    v1_.fill(0);

    Qtensor_.fill({});
}