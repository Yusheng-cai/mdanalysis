#include "GroFile.h"

void GroFile::ReadLines()
{
    ifs_.open(filename_);

    ASSERT((ifs_.is_open()), "The file with name " << filename_ << " is not opened.");

    std::string sentence;
    // The first line is comment
    std::getline(ifs_, sentence);

    // The second line is number of atoms
    std::getline(ifs_, sentence);
    num_atoms_ = StringTools::StringToType<int>(sentence); 

    // Read the main meet of the .gro file
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

    // removes the last element of the file because that is the box dimensions
    lines_.pop_back();

    // close the file as we no longer need it 
    ifs_.close();
}

void GroFile::Open(std::string filename)
{
    filename_ = filename;

    ReadLines();
    ParseFile();

    if(atomsinfo_.size() != 0)
    {
        empty_ = false;
    }
}

void GroFile::ParseFile()
{ 
    atomsinfo_.clear();
    std::set<int> ResidueSet;

    // minimum atom number
    long long minAtomNumber_ = 10000000;

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
        // residue Name is 5 characters long
        residueName = sentence.substr(5,5);
        // atomName is 5 characters long
        atomName    = sentence.substr(10,5);
        // atom Number is 5 characters long
        std::string atomNumberstr  = sentence.substr(15, 5);

        // remove the blank spaces from the string of atomName and ResidueName
        StringTools::RemoveBlankInString(residueName);
        StringTools::RemoveBlankInString(atomName);
        StringTools::RemoveBlankInString(atomName);
        StringTools::RemoveBlankInString(atomNumberstr);

        // convert the Numbers to int from string
        residueNumber = StringTools::StringToType<int>(residueNumberstr);
        atomNumber  = StringTools::StringToType<int>(atomNumberstr);

        // instantiate the Atom object and append it to atomsinfo 
        xdr::Atom a = {residueNumber, residueName, atomName, atomNumber};
        atomsinfo_.push_back(std::move(a));
        ResidueSet.insert(residueNumber);
        AtomTypes_.insert(atomName);
        ResidueNames_.insert(residueName);

        if (atomNumber < minAtomNumber_)
        {
            minAtomNumber_ = atomNumber;
        }
    }

    CorrectMinAtomNumber(minAtomNumber_);

    numResidues_ = ResidueSet.size();
    numUniqueResidues_ = ResidueNames_.size();

    // Correct for the minimum residue number
    CorrectMinResidueNumber(ResidueSet);
    constructResidues();
}

void GroFile::CorrectMinAtomNumber(int minNum)
{
    // correct for atom Number, we make it such that it always starts with 1 
    if (minNum != 1)
    {
        int diff = 1 - minNum;

        for (int i=0;i<atomsinfo_.size();i++)
        {
            auto atom = atomsinfo_[i];
            atom.atomNumber_ += diff;
        }
    }
}

void GroFile::CorrectMinResidueNumber(std::set<int>& ResidueSet)
{
    // Check for the minimum of the residue number if residue Set is not empty
    if ( ! ResidueSet.empty())
    {
        int min = *(ResidueSet.begin());

        if (min != 1)
        {
            int diff = 1-min;

            for (int i=0;i<atomsinfo_.size();i++)
            {
                auto atom = atomsinfo_[i];
                atom.residueNumber_ += diff;
            }
        }
    }
}

void GroFile::constructResidues()
{
    ResidueGroup_.resize(numResidues_);

    // construct the Residue Group
    for (int i=0;i<atomsinfo_.size();i++)
    {
        auto& A = atomsinfo_[i];
        int resNum = A.residueNumber_;

        ResidueGroup_[resNum].atoms_.push_back(A);
    }


    for (int i=0;i<ResidueGroup_.size();i++)
    {
        for (int j=1;j<ResidueGroup_[i].atoms_.size();j++)
        {
            ASSERT((ResidueGroup_[i].atoms_[j].atomNumber_ > ResidueGroup_[i].atoms_[j-1].atomNumber_), "The atom Number are not sorted.");
        }
    }
}