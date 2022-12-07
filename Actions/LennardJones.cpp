#include "LennardJones.hpp"

LennardJones::LennardJones(Real sigma, Real epsilon)
:sigma_(sigma), epsilon_(epsilon)
{
    prefactor_ = 4.0 * epsilon_;

}

LennardJones::Real LennardJones::calculate(const Real& r){
    Real s_r = sigma_ / r;
    Real s_r_6 = std::pow(s_r, 6);
    Real s_r_12= std::pow(s_r_6, 2);

    return prefactor_ * (s_r_12 - s_r_6);
}

LennardJones::Real3 LennardJones::calculate_force(const Real3& dist, const Real& r){
    Real s_r = sigma_ / r;
    Real s_r_6 = std::pow(s_r, 6);
    Real s_r_12= std::pow(s_r_6, 2);

    Real dUdr = 48 * epsilon_/r  * (0.5 * s_r_6 - s_r_12);
    Real3 drddist = dist / r;

    Real3 force = drddist * (-dUdr);

    return force;
}