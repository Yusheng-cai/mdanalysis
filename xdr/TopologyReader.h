#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "xdr/MoleculeStructs.h"

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>

// This is the topology reader
// Now we only really process the [ atoms ] directives which has all the information about the masses of the atoms

class TopologyReader
{ 
    public:
        using Real = CommonTypes::Real;

        TopologyReader(){};

        // Parse the atoms directive with the following assumptions
        // atomid atomtype resnr residuename atomname cgnr charge mass
        void Parse(std::string& name);

        // print the content
        void print();

        // make map from residue name to atom type
        void MapResnameToAtomType();

        // make the vector from indices to atomtypes
        void MapIndicesToAtom();

        std::vector<Molecule::AtomType>& getAtomTypes() {return atomtypeIndices_;}
        Molecule::atom& getAtomByIndex(int i) {return atoms_[i];}
        Molecule::residue& getResidueByIndex(int i) {return residues_[i];}

        enum TopIdx
        {
            atomnumber = 0,
            atomtype = 1,
            resname = 3,
            atomName = 4,
            charge = 6,
            mass = 7
        };
    private:
        std::string topName_;

        std::vector<std::string> resnames_;

        std::vector<Molecule::AtomType> atomtypeIndices_;

        // Each line in [ atoms ] directive must have 8 entries
        int linenum_ = 8;

        std::vector<Molecule::AtomType> atomtypes_;
        std::map<std::string, int> MapResnameToNumberResidues_;
        std::map<std::string, std::vector<Molecule::AtomType>> MapResnameToAtomType_;
        std::vector<Molecule::atom> atoms_;
        std::vector<Molecule::residue> residues_;
};