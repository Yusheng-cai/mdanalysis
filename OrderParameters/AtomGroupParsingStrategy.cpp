#include "AtomGroupParsingStrategy.h"

namespace AtomGroupParsingRegistry
{
    registry_<AtomIndexParsing> registerAtomIndex("atom_index");
    registry_<ResidueNumberParsing> registerResidueParsing("resnum");
    registry_<AtomIndexFile> registerAtomIndexFile("atom_index_file");
    registry_<AtomTypeParsing> registerAtomType("atomtype");
    registry_<ResidueNameParsing> registerResidueName("resname");
    registry_<IndexFileParsing> registerIndexFile("indexfile");
}

void AtomGroupParsingStrategy::SortAndCheckNoDuplicate(std::vector<int>& indices)
{
    // sort the indices vector
    std::sort(indices.begin(), indices.end());
    ASSERT((std::adjacent_find(indices.begin(), indices.end()) == indices.end()), "There is duplicate indices in the atom list provided."); 
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
            //ASSERT((index_str_[i].size() == 1), "Since no '-' was provided, the passed in value has to be a constant index.");

            int index = StringTools::StringToType<int>(index_str_[i]); 
            indices.push_back(index - 1);
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

                if(! grofile_.isEmpty())
                {
                    ASSERT((begin_index > 0 && end_index <= grofile_.getNumAtoms()), "Either the begin index is less than 1 or the end_index is larger than allowed. The total atoms = " << grofile_.getNumAtoms());
                }

                skip_       = StringTools::StringToType<int>(skip);
            }
            else 
            {
                std::string begin_index_str = index_str_[i].substr(0, found_dash);
                std::string end_index_str   = index_str_[i].substr(found_dash + 1);

                begin_index = StringTools::StringToType<int>(begin_index_str);
                end_index   = StringTools::StringToType<int>(end_index_str);
                if(! grofile_.isEmpty())
                {
                    ASSERT((begin_index > 0 && end_index <= grofile_.getNumAtoms()), "Either the begin index is less than 1 or the end_index is larger than allowed. The total atoms = " << grofile_.getNumAtoms());
                }

                skip_ = 1;
            }


            int j = begin_index;

            while (j <= end_index)
            {
                indices.push_back(j - 1);
                j += skip_;
            }
        }
    }

    // sort the indices vector
    SortAndCheckNoDuplicate(indices); 
}

void AtomIndexFile::Parse(std::vector<int>& indices){
    // obtain the filename
    fileName_ = selection_str_[1];

    std::ifstream ifs;
    ifs.open(fileName_);
    ASSERT((ifs.is_open()), "The file with name " << selection_str_[1] << " is not opened.");

    indices.clear();

    std::stringstream ss;

    std::string sentence;
    while(std::getline(ifs, sentence))
    {
        ss.str(sentence);
        int index;
        while (ss >> index)
        {
            indices.push_back(index-1);
        }
        ss.clear();
    }
    ifs.close(); 

    SortAndCheckNoDuplicate(indices);
}

void ResidueNumberParsing::Parse(std::vector<int>& indices)
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
                    indices.push_back(grofile_.getAtomNumber(i) - 1);
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

            while(i <= end_index)
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
                    ResidueIndices_.insert(grofile_.getResidueNumber(j)-1);
                }
            }
        }
    }

    // sort the indices vector
    SortAndCheckNoDuplicate(indices); 
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
    SortAndCheckNoDuplicate(indices); 
}

void ResidueNameParsing::Parse(std::vector<int>& indices)
{
    ASSERT((grofile_.isEmpty() == false), "The Atom Type Parsing method requires input of a gro file. However, either the gro file \
    is empty or no gro file was provided.");
    const auto& ResidueNameSet = grofile_.getResidueNames();

    for (int i =1; i<selection_str_.size();i++)
    {
        if(selection_str_[i] != ",")
        {
            index_str_.push_back(selection_str_[i]);
        }
    }

    std::vector<std::string> ResidueNames;

    for (int i=0;i<index_str_.size();i++)
    {
        auto it = ResidueNameSet.find(index_str_[i]);
        ASSERT((it != ResidueNameSet.end()), "The residue type " << index_str_[i] << " is not found"); 

        ResidueNames.push_back(index_str_[i]);
    }
 
    for (int i=0;i<grofile_.getNumAtoms();i++)
    {
        std::string resname = grofile_.getResidueName(i);

        bool found = ! (std::find(ResidueNames.begin(), ResidueNames.end(), resname) == ResidueNames.end());

        if (found)
        {
            // atom number is 1 based, while indices here are 0 based
            indices.push_back(grofile_.getAtomNumber(i) - 1);

            // insert into residue number, residue number is also 1 based
            ResidueIndices_.insert(grofile_.getResidueNumber(i) - 1);
        }
    }

    // sort the indices vector
    SortAndCheckNoDuplicate(indices);
}

IndexFileParsing::IndexFileParsing(AtomGroupParsingInput& input)
:AtomGroupParsingStrategy(input)
{
    ASSERT(( selection_str_.size() == 2), "Only 1 file name can be provided at this point.");

    // obtain the filename
    fileName_ = selection_str_[1];

    ifs_.open(fileName_);
    ASSERT((isOpen()), "The file with name " << selection_str_[1] << " is not opened.");

    std::string sentence;
    while(std::getline(ifs_, sentence))
    {
        int found = sentence.find_first_of(comment_symbol);
        if (found == std::string::npos)
        {
            ss_.str(sentence);
            int index;
            std::vector<int> index_for_sentence;

            // The first one is always the time index
            ss_ >> index;
            Frames_.push_back(index);

            while (ss_ >> index)
            {
                index_for_sentence.push_back(index);
            }

            SortAndCheckNoDuplicate(index_for_sentence);

            Fileindices_.push_back(index_for_sentence);

            ss_.clear();
        }
    }
    ifs_.close(); 

    totalFrames_ = Fileindices_.size();

    frame_count = 0;
}

bool IndexFileParsing::isOpen()
{
    return ifs_.is_open();
}

void IndexFileParsing::Parse(std::vector<int>& indices)
{
    indices.clear();

    indices.insert(indices.end(),Fileindices_[0].begin(), Fileindices_[0].end());

    //frame_count ++;
}

void IndexFileParsing::update(std::vector<int>& indices)
{
    ASSERT((frame_count < totalFrames_), "Out of range frame for Index File " << selection_str_[0]);

    indices.clear();

    indices.insert(indices.end(), Fileindices_[frame_count].begin(), Fileindices_[frame_count].end());

    frame_count++;
}
