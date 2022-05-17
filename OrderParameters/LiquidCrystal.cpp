#include "LiquidCrystal.h"

LiquidCrystal::LiquidCrystal(const OrderParametersInput& input)
:OrderParameters(input)
{
    input.pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    indicatorgroupname_ = headgroupname_;
    input.pack_.ReadString("indicatorgroup", ParameterPack::KeyType::Optional, indicatorgroupname_);
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);

    addAtomGroup(headgroupname_);
    addAtomGroup(tailgroupname_);

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    tailgroupsize_ = tailAG.getAtomGroupIndices().size();
    headgroupsize_ = headAG.getAtomGroupIndices().size();

    ASSERT((tailgroupsize_ == headgroupsize_), "The number of atoms passed in for headgroup is less than that of tail group."); 

    // register some outputs
    registerOutput("n", [this](void)->Real {return this->getN();});
    registerOutput("ntilde", [this](void)->Real {return this->getNtilde();});
}

void LiquidCrystal::CalculateDirector(std::vector<Real3>& uij, std::vector<Real>& norms, std::vector<Real>& indicators, Real& Ntilde, Real& N)
{
    auto& headAG       = getAtomGroup(headgroupname_).getAtoms();
    auto& tailAG       = getAtomGroup(tailgroupname_).getAtoms();
    auto& indicatorAG  = simstate_.getAtomGroup(indicatorgroupname_).getAtoms();
    auto& probe_volume = simstate_.getProbeVolume(pvname_);

    // resize and clear the input vectors 
    uij.clear();
    norms.clear();
    indicators.clear();
    uij.resize(headAG.size());
    norms.resize(headAG.size());
    indicators.resize(headAG.size());
    Ntilde = 0.0;
    N = 0.0;

    // start the loop
    #pragma omp parallel
    {
        Real Ntilde_local = 0.0;
        Real N_local = 0.0;
        #pragma omp for
        for (int i=0;i<headAG.size();i++)
        {
            auto pv_output = probe_volume.calculate(indicatorAG[i].position);
            indicators[i]  = pv_output.htilde_x_;

            Real3 tailatompos = tailAG[i].position;
            Real3 headatompos = headAG[i].position;

            Real3 localdirector = {};
            Real sq_dist = 0.0;
            Real norm;
            simbox_.calculateDistance(headatompos, tailatompos, localdirector, sq_dist);

            // normalize the local director
            norm = std::sqrt(sq_dist);
            localdirector = LinAlg3x3::vec_mult(1.0/norm, localdirector);

            uij[i] = localdirector;
            norms[i] = norm;

            Ntilde_local += pv_output.htilde_x_;
            N_local += pv_output.hx_;
        }

        #pragma omp critical
        {
            Ntilde += Ntilde_local;
            N      += N_local;
        }
    }
}

void LiquidCrystal::CalculateQtilde(const std::vector<Real3>& uij, Matrix& Qtensor, const std::vector<Real>& indicators, Real Ntilde)
{
    // zero out qtensor
    LinAlg3x3::matrix_zero_inplace(Qtensor);
    int numatoms = getAtomGroup(tailgroupname_).getAtoms().size();

    #pragma omp parallel
    {
        Matrix Qtensor_local;
        Qtensor_local.fill({});

        #pragma omp for
        for (int i=0;i<numatoms;i++)
        {
            auto localdirector = uij[i];

            // calculate the atomic Qtensor
            Matrix Q = LinAlg3x3::vec_dyadic(localdirector, localdirector);
            
            LinAlg3x3::matrix_mult_inplace(Q, 3);
            Q = LinAlg3x3::matrix_sub(Q, LinAlg3x3::matrix_Identity());
            LinAlg3x3::matrix_mult_inplace(Q, indicators[i]);

            // accumulate the local Qtensor
            LinAlg3x3::matrix_accum_inplace(Qtensor_local, Q);
        }

        #pragma omp critical
        {
            LinAlg3x3::matrix_accum_inplace(Qtensor, Qtensor_local);
        }
    }

    LinAlg3x3::matrix_mult_inplace(Qtensor, 1.0/(2.0*Ntilde));
}