#include "AtomGroupParsingStrategy.h"

namespace AtomGroupParsingRegistry
{
    registry_<AtomIndexParsing> registerAtomIndex("atom_index");
    registry_<ResidueParsing> registerResidueParsing("residue");
    registry_<AtomTypeParsing> registerAtomType("atomtype");
}

void AtomIndexParsing::Parse(std::vector<int>& indices)
{
    for (int i =1; i<selection_str_.size();i++)
    {
        if(selection_str_[i] != ",")
        {
            index_str_.push_back(selection_str_[i]);
        }
    }

    for ( int i=0;i<index_str_.size();i++)
    {
        int found_dash = index_str_[i].find_first_of("-");

        if (found_dash == std::string::npos)
        {
            ASSERT((index_str_[i].size() == 1), "Since no '-' was provided, the passed in value has to be a constant index.");

            int index = StringTools::StringToType<int>(index_str_[i]); 
            indices.push_back(index);
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

                ASSERT((begin_index > 0 && end_index <= grofile_.getNumAtoms()), "Either the begin index is less than 1 or the end_index is larger than allowed.");

                skip_       = StringTools::StringToType<int>(skip);
            }
            else 
            {
                std::string begin_index_str = index_str_[i].substr(0, found_dash);
                std::string end_index_str   = index_str_[i].substr(found_dash + 1);

                begin_index = StringTools::StringToType<int>(begin_index_str);
                end_index   = StringTools::StringToType<int>(end_index_str);
                ASSERT((begin_index > 0 && end_index <= grofile_.getNumAtoms()), "Either the begin index is less than 1 or the end_index is larger than allowed.");

                skip_ = 1;
            }


            int i = begin_index;

            while (i < end_index)
            {
                indices.push_back(i - 1);
                i += skip_;
            }
        }
    }

    // sort the indices vector
    std::sort(indices.begin(), indices.end());

    for (int i=0;i<indices.size()-1;i++)
    {
        ASSERT((indices[i] != indices[i+1]), "There is duplicate indices in the atom list provided.");
    }
}

void ResidueParsing::Parse(std::vector<int>& indices)
{
    ASSERT((grofile_.isEmpty() == false), "The residue Parsing method requires input of a gro file. However, either the gro file \
    is empty or no gro file was provided.");

    for (int i =1; i<selection_str_.size();i++)
    {
        if(selection_str_[i] != ",")
        {
            index_str_.push_back(selection_str_[i]);
        }
    }

    for ( int i=0;i<index_str_.size();i++)
    {
        int found_dash = index_str_[i].find_first_of("-");

        if (found_dash == std::string::npos)
        {
            ASSERT((index_str_[i].size() == 1), "Since no '-' was provided, the passed in value has to be a constant index.");

            int index = StringTools::StringToType<int>(index_str_[i]); 

            for (int i=0;i<grofile_.getNumAtoms();i++)
            {
                int residueNumber = grofile_.getResidueNumber(i);
                if(residueNumber == index)
                {
                    indices.push_back(i);
                }
            }
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

                skip_       = StringTools::StringToType<int>(skip);
            }
            else 
            {
                std::string begin_index_str = index_str_[i].substr(0, found_dash);
                std::string end_index_str   = index_str_[i].substr(found_dash + 1);

                begin_index = StringTools::StringToType<int>(begin_index_str);
                end_index   = StringTools::StringToType<int>(end_index_str);

                skip_ = 1;
            }

            ASSERT((begin_index > 0 && end_index <= grofile_.getNumResidues()), "Either the begin index is less than 1 \
             or the end_index is larger than the number of residues.");

            int i = begin_index;
            std::vector<int> list_;

            while(i < end_index)
            {
                list_.push_back(i);
                i+=skip_;
            }

            for (int j=0;j<grofile_.getNumResidues();j++)
            {
                int resNum = grofile_.getResidueNumber(j);
                bool found = (std::find(list_.begin(), list_.end(), resNum) != list_.end());

                if (found)
                {
                    indices.push_back(grofile_.getAtomNumber(j)-1);
                }
            }
        }
    }

    // sort the indices vector
    std::sort(indices.begin(), indices.end());

    for (int i=0;i<indices.size()-1;i++)
    {
        ASSERT((indices[i] != indices[i+1]), "There is duplicate indices in the atom list provided.");
    }
}

void AtomTypeParsing::Parse(std::vector<int>& indices)
{
    ASSERT((grofile_.isEmpty() == false), "The Atom Type Parsing method requires input of a gro file. However, either the gro file \
    is empty or no gro file was provided.");
    const auto& AtomTypeSet = grofile_.getAtomTypes();

    for (int i =1; i<selection_str_.size();i++)
    {
        if(selection_str_[i] != ",")
        {
            index_str_.push_back(selection_str_[i]);
        }
    }

    for (auto it = AtomTypeSet.begin();it != AtomTypeSet.end();it++)
    {
        std::cout << *it << std::endl;
    }

    std::vector<std::string> AtomTypeNames;
    for (int i=0;i<index_str_.size();i++)
    {
        auto it = AtomTypeSet.find(index_str_[i]);
        ASSERT((it != AtomTypeSet.end()), "The atom type " << index_str_[i] << " is not found"); 

        AtomTypeNames.push_back(index_str_[i]);
    }

    for (int i=0; i< grofile_.getNumAtoms();i++)
    {
        std::string atomName = grofile_.getAtomName(i);

        bool found = (std::find(AtomTypeNames.begin(), AtomTypeNames.end(), atomName) != AtomTypeNames.end());

        if (found)
        {
            indices.push_back(grofile_.getAtomNumber(i) - 1);
        }
    }

    // sort the indices vector
    std::sort(indices.begin(), indices.end());

    for (int i=0;i<indices.size()-1;i++)
    {
        ASSERT((indices[i] != indices[i+1]), "There is duplicate indices in the atom list provided.");
    }
}