#include "AtomGroup.h"

AtomGroup::AtomGroup(const ParameterPack& pack)
{
    pack.ReadVectorString("selection", ParameterPack::KeyType::Required, selection_str_);
    pack.ReadString("name", ParameterPack::KeyType::Required, name_);
 
    // the first item of selection describes the selection method --> i.e. atom_index
    for (int i =1; i<selection_str_.size();i++)
    {
        if(selection_str_[i] != ",")
        {
            index_str_.push_back(selection_str_[i]);
        }
    }

    ParseSelectionString();
}

void AtomGroup::ParseSelectionString()
{
    for ( int i=0;i<index_str_.size();i++)
    {
        int found_dash = index_str_[i].find_first_of("-");

        if (found_dash == std::string::npos)
        {
            ASSERT((index_str_[i].size() == 1), "Since no '-' was provided, the passed in value has to be a constant index.");

            int index = StringTools::StringToType<int>(index_str_[i]); 
            AtomGroupIndices_.push_back(index);
        }
        else
        {
            int found_colon = index_str_[i].find_first_of(":");
            int begin_index;
            int end_index;
            int skip_;

            if (found_colon != std::string::npos)
            {
                std::string end_index_str   = index_str_[i].substr(found_dash+1,found_colon - found_dash - 1);
                std::string begin_index_str = index_str_[i].substr(0,found_dash);
                std::string skip            = index_str_[i].substr(found_colon + 1);

                begin_index = StringTools::StringToType<int>(begin_index_str);
                end_index   = StringTools::StringToType<int>(end_index_str);
                begin_index-= 1;
                end_index  -= 1;

                skip_       = StringTools::StringToType<int>(skip);
            }
            else 
            {
                std::string begin_index_str = index_str_[i].substr(0, found_dash);
                std::string end_index_str   = index_str_[i].substr(found_dash + 1);

                begin_index = StringTools::StringToType<int>(begin_index_str);
                end_index   = StringTools::StringToType<int>(end_index_str);

                begin_index -= 1;
                end_index   -= 1;
                skip_ = 1;
            }

            int i = begin_index;

            while (i < end_index)
            {
                AtomGroupIndices_.push_back(i);
                i += skip_;
            }
        }
    }

    // sort the indices vector
    std::sort(AtomGroupIndices_.begin(), AtomGroupIndices_.end());

    for (int i=0;i<AtomGroupIndices_.size()-1;i++)
    {
        if (AtomGroupIndices_[i] == AtomGroupIndices_[i+1])
        {
            ASSERT((true == false), "There is duplicate in the atom indices passed in. Please check your atom indices again.");
        }
    }

    for (int i=0;i<AtomGroupIndices_.size();i++)
    {
        AtomGroupIndicesToGlobalIndices_.insert(std::make_pair(i, AtomGroupIndices_[i]));    
    }


    num_atoms_ = AtomGroupIndices_.size();
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

