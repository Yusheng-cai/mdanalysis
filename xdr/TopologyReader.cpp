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
        if (! sentence.empty())
        {
            std::stringstream ss;
            ss.str(sentence);

            std::string word;
            std::vector<std::string> words;

            while (ss >> word)
            {
                words.push_back(word);
            }

            if (words[0] != ";")
            {
                contents.push_back(sentence);
            }
        }
    }

    std::stringstream ss;
    int start = 0;

    for (int i=0;i<contents.size();i++)
    {
        // see if we can access [ atoms ]
        if (contents[i][0] == '[')
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
                start = i;
                break;
            }
        }
    }

    for (int i=start+1;i<contents.size();i++)
    {
        if (contents[i][0] == '[')
        {
            break;
        }
        else
        {
            ss.clear();
            ss.str(contents[i]);

            std::string word;
            std::vector<std::string> words;

            while (ss >> word)
            {
                words.push_back(word);
            }

            ASSERT((words.size() == linenum_), "The size of the line in [ atoms ] directive is " << words.size() << " while it should be " << linenum_);

            // topology file and gro file are usually needed at the same time
            AtomType a;
            a.atomName_ = words[TopIdx::atomName];
            a.charge_   = StringTools::StringToType<Real>(words[TopIdx::charge]);
            a.mass_     = StringTools::StringToType<Real>(words[TopIdx::mass]);
            a.resname_  = words[TopIdx::resname];
            a.type_     = words[TopIdx::atomtype];

            atomtypes_.push_back(a);
        }
    }        

    MakeAtomNameToMassMap();
    MakeAtomNameToTypeMap();
}

void TopologyReader::MakeAtomNameToMassMap()
{
    AtomNameToTypeMap_.clear(); 
    for (int i=0;i<atomtypes_.size();i++)
    {
        AtomNameToMassMap_.insert(std::make_pair(atomtypes_[i].atomName_, atomtypes_[i].mass_));
    }
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

TopologyReader::Real TopologyReader::getMassFromAtomName(const std::string& atomname)
{
    auto it = AtomNameToMassMap_.find(atomname);

    ASSERT((it != AtomNameToMassMap_.end()), "The atomname " << atomname << " does not exist in topology.");

    return it -> second;
}