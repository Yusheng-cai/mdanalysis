#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "MoleculeStructs.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <map>

class GroFile
{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        GroFile(){};
        ~GroFile(){};

        void Open(std::string Name);
        void ParseFile();
        void ReadLines();
        void CorrectMinResidueNumber(std::set<int>& ResidueSet);
        void constructResidues();

        // getters
        int getResidueNumber(int i)const {return atomsinfo_[i].resnum_;}
        std::string getResidueName(int i) const{return atomsinfo_[i].resname_;}
        std::string getAtomName(int i)const {return atomsinfo_[i].atomName_;}
        int getAtomNumber(int i)const{return atomsinfo_[i].atomNumber_;}
        int getNumAtoms() const {return num_atoms_;}
        int getNumResidues() const {return numResidues_;}
        int getNumUniqueResidues() const {return numUniqueResidues_;}
        const std::vector<Real3> getPosition() const {return position_;}
        const std::set<std::string> getAtomTypes() const {return AtomTypes_;}
        const std::set<std::string> getResidueNames() const {return ResidueNames_;}
        const std::vector<Molecule::residue>& getResidues() const {return ResidueGroup_;}

        // booleans that tells others whether or not the GroFile is read
        bool isEmpty() const {return empty_;}

    private:
        std::string filename_;
        std::string Cformat_ = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f";

        std::vector<Molecule::atom> atomsinfo_;

        std::ifstream ifs_;

        std::vector<std::string> lines_;

        // keeps track of the unique atomTypes in the gro file
        std::set<std::string> AtomTypes_;
        std::set<std::string> ResidueNames_;

        int numlines_;

        bool empty_=true;

        int num_atoms_;

        // Number of total residues 
        int numResidues_;

        // Number of unique residues 
        int numUniqueResidues_;

        // residues
        std::vector<Molecule::residue> ResidueGroup_;

        // residue map index to name
        std::map<int, std::string> MapIndexToResidueName_;

        // residue name map to index
        std::map<std::string, int> MapResidueNameToIndex_;

        // read the positions in gro file
        std::vector<Real3> position_;
};