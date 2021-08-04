#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

struct Atom
{
    int residueNumber_;
    std::string residueName_;
    std::string atomName_;
    int atomNumber_;
};

class GroFile
{
    public:
        GroFile(){};
        ~GroFile(){};

        void Open(std::string Name);
        void ParseFile();
        void ReadLines();

        // getters
        int getResidueNumber(int i)const {return atomsinfo_[i].residueNumber_;}
        std::string getResidueName(int i) const{return atomsinfo_[i].residueName_;}
        std::string getAtomName(int i)const {return atomsinfo_[i].atomName_;}
        int getAtomNumber(int i)const{return atomsinfo_[i].atomNumber_;}
        int getNumAtoms() const {return num_atoms_;}
        int getNumResidues() const {return numResidues_;}
        int getNumUniqueResidues() const {return numUniqueResidues_;}
        const std::set<std::string> getAtomTypes() const {return AtomTypes_;}
        const std::set<std::string> getResidueNames() const {return ResidueNames_;}

        // booleans that tells others whether or not the GroFile is read
        bool isEmpty() const {return empty_;}

    private:
        std::string filename_;
        std::string Cformat_ = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f";

        std::vector<Atom> atomsinfo_;

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
};