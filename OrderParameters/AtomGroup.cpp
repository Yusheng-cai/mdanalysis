#include "AtomGroup.h"

AtomGroup::AtomGroup(const AtomGroupInput& input)
:grofile_(input.grofile_)
{
    input.pack_.ReadVectorString("selection", ParameterPack::KeyType::Required, selection_str_);
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
 
    // the first item of selection describes the selection method --> i.e. atom_index
    AtomGroupParsingInput parsingInput = { grofile_, selection_str_ };
    strategy_ = stratptr(AtomGroupParsingRegistry::Factory::instance().create(selection_str_[0],parsingInput));
    
    // This is ensured to be sorted by AtomGroupStrategy
    strategy_->Parse(AtomGroupGlobalIndices_);
    numAtomGroupatoms_ = AtomGroupGlobalIndices_.size();

    for (int i=0;i<AtomGroupGlobalIndices_.size();i++)
    {
        AtomGroupIndicesToGlobalIndices_.insert(std::make_pair(i, AtomGroupGlobalIndices_[i]));    
        GlobalIndicesToAtomGroupIndices_.insert(std::make_pair(AtomGroupGlobalIndices_[i],i));
    }

    atoms_.resize(numAtomGroupatoms_);
}

void AtomGroup::update(const VectorReal3& total_atoms)
{
    strategy_ -> update(AtomGroupGlobalIndices_);
    numAtomGroupatoms_ = AtomGroupGlobalIndices_.size();
    atoms_.clear();
    atoms_.resize(numAtomGroupatoms_);

    for (int i=0;i < numAtomGroupatoms_;i++)
    {
        OP::Atom p;
        p.position = total_atoms[AtomGroupGlobalIndices_[i]];
        p.index = AtomGroupGlobalIndices_[i];
        atoms_[i] = p;
    }
}

int AtomGroup::AtomGroupIndices2GlobalIndices(int atomgroupIndices) const
{
    auto it = AtomGroupIndicesToGlobalIndices_.find(atomgroupIndices);

    ASSERT((it != AtomGroupIndicesToGlobalIndices_.end()), "The atomgroup indices " << atomgroupIndices << " is not found.");

    return it ->second;
}

int AtomGroup::GlobalIndices2AtomGroupIndices(int globalIndices) const 
{
    auto it = GlobalIndicesToAtomGroupIndices_.find(globalIndices);

    ASSERT((it != GlobalIndicesToAtomGroupIndices_.end()), "The global indices " << globalIndices << " is not found.");

    return it -> second;
}