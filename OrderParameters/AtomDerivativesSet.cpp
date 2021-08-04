#include "AtomDerivativesSet.h"

void AtomDerivativesSet::insert(int index, const Real3& derivatives)
{
    auto it = Index_Map_.find(index);

    ASSERT((it == Index_Map_.end()), "The index " << index << " is already in the map.");

    AtomDerivatives ad;

    ad.derivatives = derivatives;
    ad.index = index;

    AtomDerivatives_.push_back(ad);

    Index_Map_.insert(std::make_pair(index, true));
}

void AtomDerivativesSet::clear()
{
    AtomDerivatives_.clear();
    Index_Map_.clear();
}

void AtomDerivativesSet::insertOMP(int index, const Real3& derivatives)
{ 
    int numthreads = OpenMP::get_num_threads();

    if (! OpenMP::in_parallel())
    {
        ASSERT((numthreads == 1), "You are trying to insert into AtomDerivativesSet with OMP but it is neither in a parallel region, nor is number of \
        threads = 1.");

        // if OMP threads = 1, then just perform norming insertion
        insert(index, derivatives);
    }
    else
    {
        // access the atom derivatives by OMP id
        auto& atomderivative = AtomDerivativesBuffer_.access_buffer_by_id();
        
        AtomDerivatives ad;

        ad.derivatives = derivatives;
        ad.index = index;

        atomderivative.push_back(ad);
    }
}

void AtomDerivativesSet::CombineAndClearOMPBuffer()
{
    int size = AtomDerivatives_.size();
    for (auto it = AtomDerivativesBuffer_.beginworker(); it != AtomDerivativesBuffer_.endworker();it++)
    {
        size += it -> size();
    }


    AtomDerivatives_.reserve(size);

    for (auto it = AtomDerivativesBuffer_.beginworker();it!=AtomDerivativesBuffer_.endworker();it++)
    {
        AtomDerivatives_.insert(AtomDerivatives_.end(), it->begin(), it ->end());
    }

    // clear the buffer
    AtomDerivativesBuffer_.clearBuffer();
}

void AtomDerivativesSet::setMasterObject()
{
    AtomDerivativesBuffer_.set_master_object(AtomDerivatives_);
}