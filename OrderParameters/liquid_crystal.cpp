#include "liquid_crystal.h"

liquid_crystal::liquid_crystal(const OrderParametersInput& input)
:OrderParameters(input)
{
    input.pack_.ReadString("tailgroup", ParameterPack::KeyType::Required, tailgroupname_);
    input.pack_.ReadString("headgroup", ParameterPack::KeyType::Required, headgroupname_);

    addAtomGroup(headgroupname_);
    addAtomGroup(tailgroupname_);

    auto& tailAG = getAtomGroup(tailgroupname_);
    auto& headAG = getAtomGroup(headgroupname_);

    tailgroupsize_ = tailAG.getAtomGroupIndices().size();
    headgroupsize_ = headAG.getAtomGroupIndices().size();

    ASSERT((tailgroupsize_ == headgroupsize_), "The number of atoms passed in for headgroup is less than that of tail group."); 
}

void liquid_crystal::getUij()
{
    uij_.clear();
    norms_.clear();

    auto& headAG = getAtomGroup(headgroupname_);
    auto& tailAG = getAtomGroup(tailgroupname_);

    auto& headatoms_ = headAG.getAtoms();
    auto& tailatoms_ = tailAG.getAtoms();

    ASSERT((tailatoms_.size() == headatoms_.size()), "The size of the tail atom is " << tailatoms_.size() << " while the size of the head atom is \
    " << headatoms_.size());

    uij_.resize(headatoms_.size());
    norms_.resize(headatoms_.size());

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
            Real norm;
            simbox_.calculateDistance(headatompos_, tailatompos_, localdirector, sq_dist);

            // normalize the local director
            norm = std::sqrt(sq_dist);
            localdirector = Qtensor::vec_mult(1.0/norm, localdirector);

            uij_[i] = localdirector;
            norms_[i] = norm;
        }
    }
}

void liquid_crystal::calcQtensor()
{
    Qtensor_.fill({});

    auto& headAG = getAtomGroup(headgroupname_);
    auto& tailAG = getAtomGroup(tailgroupname_);

    auto& headatoms_ = headAG.getAtoms();
    auto& tailatoms_ = tailAG.getAtoms();

    #pragma omp parallel
    {
        Matrix Qtensor_local;
        Qtensor_local.fill({});

        #pragma omp for
        for (int i=0;i<tailatoms_.size();i++)
        {
            auto localdirector = uij_[i];

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
}