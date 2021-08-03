#include "AtomGroup.h"

AtomGroup::AtomGroup(const AtomGroupInput& input)
:grofile_(input.grofile_)
{
    input.pack_.ReadVectorString("selection", ParameterPack::KeyType::Required, selection_str_);
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
 
    // the first item of selection describes the selection method --> i.e. atom_index
    AtomGroupParsingInput parsingInput = { grofile_, selection_str_ };
    strategy_ = AtomGroupParsingRegistry::Factory::instance().create(selection_str_[0],parsingInput);

    strategy_->Parse(AtomGroupIndices_);
    num_atoms_ = AtomGroupIndices_.size();

    for (int i=0;i<AtomGroupIndices_.size();i++)
    {
        AtomGroupIndicesToGlobalIndices_.insert(std::make_pair(i, AtomGroupIndices_[i]));    
    }
}

void AtomGroup::update(const VectorReal3& total_atoms)
{
    atoms_positions_buffer_.clearBuffer();
    atom_positions_.clear();

    // important! address of a vector changes throughout
    atoms_positions_buffer_.set_master_object(atom_positions_);

    #pragma omp parallel
    {
        auto& buffer_ = atoms_positions_buffer_.access_buffer_by_id();
        buffer_.resize(0);

        #pragma omp for
        for (int i=0;i < num_atoms_;i++)
        {
            auto& a = total_atoms[AtomGroupIndices_[i]];
            buffer_.push_back(a);
        }
    }

    atom_positions_.reserve(num_atoms_);

    for(auto it= atoms_positions_buffer_.beginworker();it != atoms_positions_buffer_.endworker();it++)
    {
        atom_positions_.insert(atom_positions_.end(), it->begin(), it->end());
    }
}

int AtomGroup::AtomGroupIndices2GlobalIndices(int atomgroupIndices) const
{
    auto it = AtomGroupIndicesToGlobalIndices_.find(atomgroupIndices);

    ASSERT((it != AtomGroupIndicesToGlobalIndices_.end()), "The indices of number " << atomgroupIndices << " is not found.");

    return it ->second;
}

