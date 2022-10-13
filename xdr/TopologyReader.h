#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "xdr/MoleculeStructs.h"
#include "tools/FileSystem.h"
#include "tools/CommonOperations.h"

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <algorithm>

class TopologyReader
{ 
    public:
        using Real = CommonTypes::Real;

        TopologyReader(){};

        // Parse the atoms directive with the following assumptions
        // atomid atomtype resnr residuename atomname cgnr charge mass
        void Parse(std::string& name);

        // read the file
        void ReadFile(std::string& name, std::vector<std::string>& contents);

        // make map from residue name to atom type
        void MapResnameToAtomType();

        // make the vector from indices to atomtypes
        void MapIndicesToAtom();

        Molecule::atom& getAtomByIndex(int i) {return atoms_[i];}
        Molecule::residue& getResidueByIndex(int i) {return residues_[i];}


        enum TopologyItemIndex{
            atomnumber = 0,
            atomtype = 1,
            resname = 3,
            atomName = 4,
            charge = 6,
            mass = 7
        };
    private:
        std::string topName_;

        std::vector<std::string> unique_resnames_;
        std::vector<std::string> resnames_;

        std::vector<std::string> comment_str_={";"};

        std::vector<int> AtomtypeIndices_;
        std::vector<int> MoleculeTypeIndices_;
        std::vector<int> AtomIndices_;

        std::map<std::string, Molecule::AtomType> MapTypenameToAtomType_;
        std::map<std::string, int> MapResnameToNumberResidues_;
        std::map<std::string, std::vector<std::string>> MapResnameToTypename_;
        std::map<std::string, std::vector<std::string>> MapResnameToAtomname_;

        std::vector<Molecule::atom> atoms_;
        std::vector<Molecule::residue> residues_;

        // the abs path
        std::string absolute_path_;

        // if parse has been run
        bool is_empty = true;
};