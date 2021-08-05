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
    // atoms_.clear();
    // atoms_.reserve(numAtomGroupatoms_);
    // atoms_buffer_.clearBuffer();

    // // important! address of a vector changes throughout
    // atoms_buffer_.set_master_object(atoms_);

    for (int i=0;i < numAtomGroupatoms_;i++)
    {
        OP::Atom p;
        p.position = total_atoms[AtomGroupGlobalIndices_[i]];
        p.index = AtomGroupGlobalIndices_[i];
        atoms_[i] = p;
    }

    // #pragma omp parallel
    // {
    //     auto& abuffer_ = atoms_buffer_.access_buffer_by_id();
    //     abuffer_.resize(0);

    //     #pragma omp for schedule(static)
    //     for (int i=0;i < numAtomGroupatoms_;i++)
    //     {
    //         OP::Atom p;
    //         p.position = total_atoms[AtomGroupGlobalIndices_[i]];
    //         p.index = AtomGroupGlobalIndices_[i];
    //         abuffer_.push_back(p);
    //     }
    // }

    // if (OpenMP::get_max_threads() > 1)
    // {
    //     atoms_.reserve(numAtomGroupatoms_);

    //     for(auto it= atoms_buffer_.beginworker();it != atoms_buffer_.endworker();it++)
    //     {
    //         atoms_.insert(atoms_.end(), it->begin(), it->end());
    //     }
    // }

    // // want to add an extra check just to make sure that the atoms are in order
    // if (OpenMP::get_max_threads() > 1)
    // {
    //     std::cout << "We checking order" << std::endl;
    //     #pragma omp parallel for
    //     for (int i=1;i<atoms_.size();i++)
    //     {
    //         ASSERT((atoms_[i].index >= atoms_[i-1].index), "The atoms indices are not in order!");
    //     }
    // }
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