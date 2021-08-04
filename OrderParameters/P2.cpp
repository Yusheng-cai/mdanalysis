#include "P2.h"

namespace OrderParametersRegistry
{
    registry_<P2> registerP2("p2");
}

P2::P2(const OrderParametersInput& input)
:OrderParameters(input)
{
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);

    addAtomGroup(headgroupname_);
    addAtomGroup(tailgroupname_);

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    tailgroupsize_ = tailAG.getAtomGroupIndices().size();
    headgroupsize_ = headAG.getAtomGroupIndices().size();

    ASSERT((tailgroupsize_ == headgroupsize_), "The number of atoms passed in for headgroup is less than that of tail group.");

    registerOutput("p2",[this](void)-> Real {return this->getP2();});
    registerOutput("qxx", [this](void)->Real {return this->getQxx();});
    registerOutput("qxy", [this](void)->Real {return this->getQxy();});
    registerOutput("qxz", [this](void)->Real {return this ->getQxz();});
    registerOutput("qyy", [this](void)->Real {return this->getQyy();});
    registerOutput("qyz", [this](void) ->Real {return this->getQyz();});
}

void P2::calculate()
{
    // zero the Qtensor
    Qtensor_.fill({});

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    auto& tailatoms_ = tailAG.getAtoms();
    auto& headatoms_ = headAG.getAtoms();

    #pragma omp parallel
    {
        Matrix Qtensor_local = {};
        #pragma omp for
        for (int i=0;i<tailatoms_.size();i++)
        {
            Real3 tailatompos_ = tailatoms_[i].position;
            Real3 headatompos_ = headatoms_[i].position;

            Real3 localdirector = {};
            Real sq_dist = 0.0;
            simbox_.calculateDistance(headatompos_, tailatompos_, localdirector, sq_dist);

            // normalize the local director
            localdirector = Qtensor::vec_mult(1.0/std::sqrt(sq_dist), localdirector);

            // calculate the atomic Qtensor
            Matrix Qtensor_atomic = Qtensor::vec_dyadic(localdirector, localdirector);
            
            Qtensor::matrix_mult_inplace(Qtensor_atomic, 3);
            Qtensor_atomic = Qtensor::matrix_sub(Qtensor_atomic, Qtensor::matrix_Identity());

            // accumulate the local Qtensor
            Qtensor::matrix_accum_inplace(Qtensor_local, Qtensor_atomic);
        }

        #pragma omp critical
        {
            Qtensor::matrix_accum_inplace(Qtensor_, Qtensor_local);
        }
    }

    Qtensor::matrix_mult_inplace(Qtensor_, 1.0/(2.0*tailgroupsize_));

    // find the p2 variable as well as the eigenvalue
    auto pair = Qtensor::OP_Qtensor(Qtensor_);
    P2_OP_ = pair.first;
    v1_    = pair.second;
}

void P2::update()
{
    P2_OP_ = 0;
    v1_.fill(0);

    Qtensor_.fill({});
}