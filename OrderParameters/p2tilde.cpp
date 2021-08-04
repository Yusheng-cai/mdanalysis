#include "P2tilde.h"

namespace OrderParametersRegistry
{
    registry_<P2tilde> registerP2tilde("p2tilde");
}

P2tilde::P2tilde(const OrderParametersInput& input)
:OrderParameters(input), pv_(const_cast<ProbeVolumeRegistry&>(input.pv_registry_))
{
    input.pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);

    registerOutput("p2tilde", [this](void)->Real {return this->getP2tilde();});
    registerOutput("qxx", [this](void) -> Real {return this->getQxx();});
    registerOutput("qxy", [this](void) -> Real {return this->getQxy();});
    registerOutput("qxz", [this](void) -> Real {return this -> getQxz();});
    registerOutput("qyy", [this](void)-> Real {return this -> getQyy();});
    registerOutput("qyz", [this](void) -> Real {return this -> getQyz();});
}

void P2tilde::calculate()
{
    auto& probeV = pv_.getProbeVolume(pvname_);
    auto& headAG = simstate_.getAtomGroup(headgroupname_);
    auto& tailAG = simstate_.getAtomGroup(tailgroupname_);

    auto& headatoms_ = headAG.getAtoms();
    auto& tailatoms_ = tailAG.getAtoms();

    #pragma omp parallel
    {
        auto& buffer_ = IndusIndicesbuffer_.access_buffer_by_id(); 
        Matrix Qtensor_local_ = {};
        Real localNtilde_ = 0.0;
        Real localN = 0;

        #pragma omp for
        for (int i=0; i<headatoms_.size();i++)
        {
            Real3 headpos_ = headatoms_[i].position;
            Real3 tailpos_ = tailatoms_[i].position; 

            ProbeVolumeOutput output = probeV.calculate(headpos_);

            // perform the Qtensor calculation only if htilde_x is larger than 0
            if (output.htilde_x_ > 0)
            {
                buffer_.push_back(headAG.AtomGroupIndices2GlobalIndices(i));

                Real3 local_director;
                Real sq_dist;
                simbox_.calculateDistance(headpos_, tailpos_, local_director, sq_dist);

                local_director = Qtensor::vec_mult(local_director, 1.0/std::sqrt(sq_dist));

                Matrix Qtensor_atomic = Qtensor::vec_dyadic(local_director, local_director);

                Qtensor::matrix_mult_inplace(Qtensor_atomic, 3);

                Qtensor_atomic = Qtensor::matrix_sub(Qtensor_atomic, Qtensor::matrix_Identity());
                Qtensor::matrix_mult_inplace(Qtensor_atomic, output.htilde_x_);

                Qtensor::matrix_accum_inplace(Qtensor_local_, Qtensor_atomic);

                localNtilde_ += output.htilde_x_;
                localN       += output.hx_;
            }
        }

        #pragma omp critical
        {
            Qtensor::matrix_accum_inplace(Qtensor_, Qtensor_local_);

            Ntilde_ += localNtilde_;
            N_      += localN;
        }
    }

    Qtensor::matrix_mult_inplace(Qtensor_, 1.0/(2.0*Ntilde_));
}

void P2tilde::update()
{
    auto& probeV = pv_.accessProbeVolume(pvname_);
    probeV.update();

    Qtensor_.fill({});
    p2tilde_ = 0.0;

    Ntilde_ = 0.0;
    N_ = 0.0;
}