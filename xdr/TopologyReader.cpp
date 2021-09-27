#include "TopologyReader.h"

void TopologyReader::Parse(std::string& name)
{
    atomtypes_.clear();
    
    // First read everything in the file
    std::ifstream ifs_;

    ifs_.open(name);
    ASSERT((ifs_.is_open()), "The name " << name << " is not opened.");

    std::vector<std::string> contents;
    std::string sentence;

    while (std::getline(ifs_,sentence))
    {
        // skip the empty lines
        if (! StringTools::CheckIfOnlyWhiteSpace(sentence))
        {
            std::stringstream ss;
            ss.str(sentence);

            std::string word;
            std::vector<std::string> words;

            while (ss >> word)
            {
                // append to words if word is not just empty space
                if (! StringTools::CheckIfOnlyWhiteSpace(word))
                {
                    words.push_back(word);
                }
            }

            if (words[0] != ";" && words[0] != "#include")
            {
                contents.push_back(sentence);
            }
        }
    }


    std::stringstream ss;
    std::vector<int> startIndices_;
    int start = 0;

    for (int i=0;i<contents.size();i++)
    {
        // see if we can access [ atoms ]
        // find the first one not space 
        int index = contents[i].find_first_not_of(" ");
        if (contents[i][index] == '[')
        {
            // set the stringstream
            ss.clear();
            ss.str(contents[i]);

            std::string word;
            std::vector<std::string> words;

            while (ss >> word)
            {
                words.push_back(word);
            }

            // if the word is atoms, then we start parsing
            if (words[1] == "atoms")
            {
                startIndices_.push_back(i);
            }
        }
    }
    for (int i =0;i<startIndices_.size();i++)
    {
        int start = startIndices_[i];
        for (int j=start+1;j<contents.size();j++)
        {
            int index = contents[j].find_first_not_of(" ");
            if (contents[j][index] == '[')
            {
                break;
            }
            else
            {
                ss.clear();
                ss.str(contents[j]);

                std::string word;
                std::vector<std::string> words;

                while (ss >> word)
                {
                    words.push_back(word);
                }

                ASSERT((words.size() == linenum_), "The size of the line in [ atoms ] directive is " << words.size() << " while it should be " << linenum_);

                AtomType a;
                a.atomName_ = words[TopIdx::atomName];
                a.charge_   = StringTools::StringToType<Real>(words[TopIdx::charge]);
                a.mass_     = StringTools::StringToType<Real>(words[TopIdx::mass]);
                a.resname_  = words[TopIdx::resname];
                a.type_     = words[TopIdx::atomtype];

                atomtypes_.push_back(a);
            }
        }        
    }

    MakeResnameAtomNameToTypeMap();
    MakeResnameAtomTypeToMassMap();
    MakeResnameAtomTypeToChargeMap();
}

void TopologyReader::MakeResnameAtomNameToTypeMap()
{
    ResNameAtomNameToTypeMap_.clear();

    for (int i=0;i<atomtypes_.size();i++)
    {
        std::vector<std::string> str_vec_(2);

        str_vec_[0] = atomtypes_[i].resname_;
        str_vec_[1] = atomtypes_[i].atomName_;

        auto it = ResNameAtomNameToTypeMap_.find(str_vec_);

        if (it  == ResNameAtomNameToTypeMap_.end())
        {
            ResNameAtomNameToTypeMap_.insert(std::make_pair(str_vec_, atomtypes_[i].type_));
        }
    }
}

void TopologyReader::MakeResnameAtomTypeToChargeMap()
{
    ResNameAtomTypeToChargeMap_.clear();

    for (int i=0;i<atomtypes_.size();i++)
    {
        std::vector<std::string> str_vec_(2);
        str_vec_[0] = atomtypes_[i].resname_;
        str_vec_[1] = atomtypes_[i].type_;

        auto it = ResNameAtomTypeToChargeMap_.find(str_vec_);

        if (it == ResNameAtomTypeToChargeMap_.end())
        {
            ResNameAtomTypeToChargeMap_.insert(std::make_pair(str_vec_, atomtypes_[i].charge_));
        }
    }
}

TopologyReader::Real TopologyReader::getChargeFromAtomTypeResname(const std::string& resname, const std::string& atomType)
{
    std::vector<std::string> str_vec_(2);

    str_vec_[0] = resname;
    str_vec_[1] = atomType;

    auto it = ResNameAtomTypeToChargeMap_.find(str_vec_);

    ASSERT((it != ResNameAtomTypeToChargeMap_.end()), "The atom with resname " << resname << " and type " << atomType << " does not exist in topology.");

    return it -> second;
}

void TopologyReader::MakeResnameAtomTypeToMassMap()
{
    ResNameAtomTypeToMassMap_.clear(); 

    for (int i=0;i<atomtypes_.size();i++)
    {
        std::vector<std::string> str_vec_(2);
        str_vec_[0] = atomtypes_[i].resname_;
        str_vec_[1] = atomtypes_[i].type_;
        auto it  = ResNameAtomTypeToMassMap_.find(str_vec_);

        if (it == ResNameAtomTypeToMassMap_.end())
        {
            ResNameAtomTypeToMassMap_.insert(std::make_pair(str_vec_,atomtypes_[i].mass_));
        }
    }
}

std::string TopologyReader::getAtomTypeFromAtomNameResname(const std::string& resname, const std::string& atomName)
{   
    std::vector<std::string> str_vec_(2);
    str_vec_[0] = resname;
    str_vec_[1] = atomName;

    auto it = ResNameAtomNameToTypeMap_.find(str_vec_);

    ASSERT((it != ResNameAtomNameToTypeMap_.end()), "The atomname " << atomName << " and resname " << resname << " does not exist in topology.");

    return it -> second;

}

void TopologyReader::print()
{
    std::cout << "atomName\tcharge\tmass\tresname\ttype" << std::endl;

    for(int i=0;i<atomtypes_.size();i++)
    {
        auto& a = atomtypes_[i];

        std::cout << a.atomName_ << "\t" << a.charge_ << "\t" << a.mass_ << "\t" << a.resname_ << "\t" << a.type_ << std::endl;
    }
}

TopologyReader::Real TopologyReader::getMassFromAtomTypeResname(const std::string& resname, const std::string& atomType)
{
    std::vector<std::string> str_vec_(2);
    str_vec_[0] = resname;
    str_vec_[1] = atomType;
    auto it = ResNameAtomTypeToMassMap_.find(str_vec_);

    ASSERT((it != ResNameAtomTypeToMassMap_.end()), "The atomtype " << atomType << " with resname " << resname << " does not exist in topology.");

    return it -> second;
}