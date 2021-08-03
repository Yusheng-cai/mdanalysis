#include "GroFile.h"

void GroFile::ReadNumLines()
{
    ifs_.open(filename_);

    ASSERT((ifs_.is_open()), "The file with name " << filename_ << " is not opened.");

    std::string sentence;
    std::getline(ifs_, sentence);

    std::getline(ifs_, sentence);
    num_atoms_ = StringTools::StringToType<int>(sentence); 

    while ( ! ifs_.eof() )
    {
        std::string sentence;
        std::getline(ifs_, sentence);

        if (! sentence.empty())
        {
            numlines_++;
            lines_.push_back(sentence); 
        }
    }

    // removes the last element of the file
    lines_.pop_back();

    ifs_.close();
}

void GroFile::Open(std::string filename)
{
    filename_ = filename;

    ReadNumLines();
    ParseFile();

    if(atomsinfo_.size() != 0)
    {
        empty_ = false;
    }
}

void GroFile::ParseFile()
{ 
    atomsinfo_.clear();
    long long minAtomNumber = 1000000000;
    std::set<int> ResidueSet;

    for (int i=0; i< lines_.size() ;i++)
    {
        std::string sentence = lines_[i];
        ASSERT((! sentence.empty()), "The sentence read is empty.");

        int residueNumber;
        std::string residueName;
        std::string atomName;
        int atomNumber;

        // residue Number is 5 characters long
        std::string residueNumberstr = sentence.substr(0,5);
        residueName = sentence.substr(5,5);
        atomName    = sentence.substr(10,5);

        // remove the blank spaces from the string of atomName and ResidueName
        residueName.erase(std::remove(residueName.begin(),residueName.end(),' '),residueName.end());
        atomName.erase(std::remove(atomName.begin(),atomName.end(),' '),atomName.end());

        std::string atomNumberstr  = sentence.substr(15, 5);

        residueNumber = StringTools::StringToType<int>(residueNumberstr);
        atomNumber  = StringTools::StringToType<int>(atomNumberstr);

        Atom a = {residueNumber, residueName, atomName, atomNumber};
        atomsinfo_.push_back(std::move(a));
        ResidueSet.insert(residueNumber);
        AtomTypes_.insert(atomName);

        if (atomNumber < minAtomNumber)
        {
            minAtomNumber = atomNumber;
        }
    }

    numResidues_ = ResidueSet.size();

    if (minAtomNumber != 1)
    {
        int diff = 1 - minAtomNumber;

        for (int i=0;i<atomsinfo_.size();i++)
        {
            auto atom = atomsinfo_[i];
            atom.atomNumber_ += diff;
        }
    }
}