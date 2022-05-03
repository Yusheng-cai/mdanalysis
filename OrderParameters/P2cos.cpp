#include "P2cos.h"

namespace OrderParametersRegistry
{
    registry_<P2cos> registerP2cos("p2cos");
}

P2cos::P2cos(const OrderParametersInput& input)
:LiquidCrystal(input)
{
    input.pack_.ReadArrayNumber<Real,3>("director", ParameterPack::KeyType::Required, n_);

    LinAlg3x3::normalize(n_);

    registerOutput("p2",[this](void)-> Real {return this->getP2cos();});
}

void P2cos::calculate()
{
    getUij();
    P2cos_OP_ = 0.0;

    #pragma omp parallel
    {
        Real P2cos_local_ = 0.0;
        #pragma omp for
        for (int i=0;i<tailgroupsize_;i++)
        {
            Real3 AtomDirector = uij_[i];

            Real dot_product = LinAlg3x3::vec_dot(AtomDirector, n_);

            P2cos_local_ += 1.5*std::pow(dot_product,2.0) - 0.5;
        }

        #pragma omp critical
        {
            P2cos_OP_  += P2cos_local_;
        }
    }

    P2cos_OP_ = 1.0/tailgroupsize_*P2cos_OP_;

    // clear derivatives
    clearDerivativesOutputs();

    auto& headAG = getAtomGroup(headgroupname_);
    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headatoms_ = headAG.getAtoms();
    auto& tailatoms_ = tailAG.getAtoms();

    auto& head_derivatives = accessDerivatives(headgroupname_);
    auto& tail_derivatives = accessDerivatives(tailgroupname_);
    Real N = tailgroupsize_;

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<tailgroupsize_;i++)
        {
            auto headatom = headatoms_[i];
            auto tailatom = tailatoms_[i];

            auto director = uij_[i];
            auto norm = norms_[i];

            auto derivatives = dP2cosdr(N, norm, director, n_);

            head_derivatives.insertOMP(headatom.index, derivatives.first);
            tail_derivatives.insertOMP(tailatom.index, derivatives.second);
        }
    }

    head_derivatives.CombineAndClearOMPBuffer();
    tail_derivatives.CombineAndClearOMPBuffer();
}

void P2cos::update()
{

}

std::pair<P2cos::Real3,P2cos::Real3> P2cos::dP2cosdr(Real N, Real norm, Real3& director, Real3& n)
{
    std::pair<Real3,Real3> output;

    Real dot_product = LinAlg3x3::vec_dot(director, n);
    Real factor = 3.0/(N*norm)*dot_product;
    Real3 derivative;

    for (int i=0;i<3;i++)
    {
        derivative[i] = n[i] - director[i]*dot_product;
    }

    derivative = LinAlg3x3::vec_mult(factor, derivative);

    Real3 secondD = LinAlg3x3::vec_mult(-1.0, derivative);

    output.first = derivative;
    output.second = secondD;

    return output;
}