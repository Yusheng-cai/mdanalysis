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

    MakeAtomTypeToMassMap();
    MakeAtomNameToTypeMap();
    MakeAtomTypeToChargeMap();
}

void TopologyReader::MakeAtomTypeToMassMap()
{
    AtomNameToTypeMap_.clear(); 
    for (int i=0;i<atomtypes_.size();i++)
    {
        auto it  = AtomTypeToMassMap_.find(atomtypes_[i].type_);

        if (it == AtomTypeToMassMap_.end())
        {
            AtomTypeToMassMap_.insert(std::make_pair(atomtypes_[i].type_, atomtypes_[i].mass_));
        }
    }
}

void TopologyReader::MakeAtomTypeToChargeMap()
{
    AtomTypeToChargeMap_.clear();

    for (int i=0;i<atomtypes_.size();i++)
    {
        auto it = AtomTypeToChargeMap_.find(atomtypes_[i].atomName_);

        if (it == AtomTypeToChargeMap_.end())
        {
            AtomTypeToChargeMap_.insert(std::make_pair(atomtypes_[i].atomName_, atomtypes_[i].charge_));
        }
    }
    
}

std::string TopologyReader::getAtomTypeFromAtomName(const std::string& atomName)
{   
    auto it = AtomNameToTypeMap_.find(atomName);

    ASSERT((it != AtomNameToTypeMap_.end()), "The atomname " << atomName << " does not exist in topology.");

    return it -> second;

}

void TopologyReader::MakeAtomNameToTypeMap()
{
    AtomNameToTypeMap_.clear();
    for (int i=0;i<atomtypes_.size();i++)
    {
        AtomNameToTypeMap_.insert(std::make_pair(atomtypes_[i].atomName_, atomtypes_[i].type_));
    }
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

TopologyReader::Real TopologyReader::getMassFromAtomType(const std::string& atomType)
{
    auto it = AtomTypeToMassMap_.find(atomType);

    ASSERT((it != AtomTypeToMassMap_.end()), "The atomname " << atomType << " does not exist in topology.");

    return it -> second;
}

TopologyReader::Real TopologyReader::getChargeFromAtomType(const std::string& atomType)
{
    auto it  = AtomTypeToChargeMap_.find(atomType);

    ASSERT((it != AtomTypeToChargeMap_.end()), "The atom type " << atomType << " does not exist in topology.");

    return it -> second;
}