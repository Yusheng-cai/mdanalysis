#include "MorsePotential.hpp"

MorsePotential::MorsePotential(Real De, Real re, Real alpha)
:De_(De), re_(re), alpha_(alpha)
{

}

MorsePotential::Real MorsePotential::calculate(const Real& r){
    Real diff = r - re_;

    return De_ * ( std::exp(- 2 * alpha_ * diff) - 2 * std::exp(- alpha_ * diff));
}

MorsePotential::Real3 MorsePotential::calculate_force(const Real3& dist, const Real& r){
    Real diff = r - re_;
    Real dUdr = 2.0 * alpha_ * De_ * ( std::exp(-alpha_ * diff) - std::exp(-2 * alpha_ * diff) );
    Real3 drddist = dist / r;

    Real3 force = drddist * (-dUdr);

    return force;
}
