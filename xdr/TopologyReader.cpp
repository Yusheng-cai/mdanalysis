#include "TopologyReader.h"

void TopologyReader::Parse(std::string& name)
{
    atomtypes_.clear();
    
    // First read everything in the file
    std::ifstream ifs;

    // open the file and make sure it is opened
    ifs.open(name);
    ASSERT((ifs.is_open()), "The name " << name << " is not opened.");

    // start parsing the file 
    std::vector<std::string> contents;
    std::string sentence;

    // start parsing the sentences 
    while (std::getline(ifs,sentence))
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

            // no # symbol in the first letter of the sentence 
            bool NotComment = words[0].find("#") == std::string::npos;
            if (words[0] != ";" && NotComment)
            {
                contents.push_back(sentence);
            }
        }
    }

    // contents is the content of the file --> all the lines inside the file 
    std::stringstream ss;
    std::vector<int> StartIndices;
    int moleculeIndices;
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
                StartIndices.push_back(i);
            }

            if (words[1] == "molecules")
            {
                moleculeIndices = i;
            }
        }
    }

    // read the molecules
    for (int j = moleculeIndices+1;j<contents.size();j++)
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
                if (word == ";")
                {
                    // break out of while loop
                    break;
                }
                words.push_back(word);
            }

            ASSERT((words.size() == 2), "For molecule directives, we must only have 2 entries, [molecuname, number], while we have " << word.size());

            auto it = MapResnameToNumber_.find(words[0]);

            ASSERT((it == MapResnameToNumber_.end()), "The molecule name " << words[0] << " appeared twice.");

            int num = StringTools::StringToType<int>(words[1]);

            MapResnameToNumber_.insert(std::make_pair(words[0], num));

            resnames_.push_back(words[0]);
        }
    }

    // read the starting indices
    for (int i =0;i<StartIndices.size();i++)
    {
        int start = StartIndices[i];
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
                    // sometimes we have comments at the end of line
                    // e.g. blah blah ; comment 
                    if (word == ";")
                    {
                        break;
                    }
                    words.push_back(word);
                }

                ASSERT((words.size() == linenum_), "The size of the line in [ atoms ] directive is " << words.size() << " while it should be " << linenum_);

                Molecule::AtomType a;
                a.atomName_ = words[TopIdx::atomName];
                a.charge_   = StringTools::StringToType<Real>(words[TopIdx::charge]);
                a.mass_     = StringTools::StringToType<Real>(words[TopIdx::mass]);
                a.resname_  = words[TopIdx::resname];
                a.type_     = words[TopIdx::atomtype];
                a.atomNumber_ = StringTools::StringToType<int>(words[TopIdx::atomnumber]);

                atomtypes_.push_back(a);
            }
        }        
    }

    // map index to corresponding atom type information
    MakeResidueToAtomTypeMap();
    MakeIndicesToAtomType();
}

void TopologyReader::MakeIndicesToAtomType()
{
    for (int i=0;i<resnames_.size();i++)
    {
        std::string name = resnames_[i];

        auto it = MapResnameToNumber_.find(name);

        ASSERT((it != MapResnameToNumber_.end()), "The resname " << name << " is not found.");
        int num = it -> second;

        auto atit = MapResidueToAtomType_.find(name);
        ASSERT((atit != MapResidueToAtomType_.end()), "The resname " << name << " is not found.");

        auto atomtypevec = atit -> second;

        for (int j=0;j<num;j++)
        {
            atomtypeIndices_.insert(atomtypeIndices_.end(), atomtypevec.begin(), atomtypevec.end());
        }
    }
}

TopologyReader::Real TopologyReader::getChargeFromIndex(int index)
{
    ASSERT((index < atomtypeIndices_.size()), "Index " << index << " out of range.");
    return atomtypeIndices_[index].charge_;
}

TopologyReader::Real TopologyReader::getMassFromIndex(int index)
{
    ASSERT((index < atomtypeIndices_.size()), "Index " << index << " out of range.");
    return atomtypeIndices_[index].mass_;
}

std::string TopologyReader::getAtomTypeFromIndex(int index)
{
    ASSERT((index < atomtypeIndices_.size()), "Index " << index << " out of range.");
    return atomtypeIndices_[index].type_;
}

void TopologyReader::MakeResidueToAtomTypeMap()
{
    for (int i=0;i<atomtypes_.size();i++)
    {
        std::string resname = atomtypes_[i].resname_;

        auto it = MapResidueToAtomType_.find(resname);

        if (it != MapResidueToAtomType_.end())
        {
            it -> second.push_back(atomtypes_[i]);
        }
        else
        {
            std::vector<Molecule::AtomType> temp = {atomtypes_[i]};
            MapResidueToAtomType_.insert(std::make_pair(resname, temp));
        }
    }

    // sort the atomtypes in the map by their atomnumber 
    for (auto it = MapResidueToAtomType_.begin(); it != MapResidueToAtomType_.end(); it ++)
    {
        auto types = it -> second;
        std::sort(types.begin(), types.end(), [](const Molecule::AtomType& a, const Molecule::AtomType& b)\
        {return a.atomNumber_> b.atomNumber_;});
        it -> second = types;
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
