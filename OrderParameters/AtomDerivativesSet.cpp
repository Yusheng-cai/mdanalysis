#include "AtomDerivativesSet.h"

void AtomDerivativesSet::insert(int index, const Real3& derivatives)
{
    auto it = Index_Map_.find(index);

    ASSERT((it == Index_Map_.end()), "The index " << index << " is already in the map.");

    AtomDerivatives ad;

    ad.derivatives = derivatives;
    ad.index = index;

    AtomDerivatives_.push_back(ad);
}

void AtomDerivativesSet::clear()
{
    AtomDerivatives_.clear();
    Index_Map_.clear();
}