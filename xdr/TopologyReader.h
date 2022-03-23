#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

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

        // make map from atomName to mass
        void MakeResnameAtomTypeToMassMap();

        // make map from resname and atom name to charge --> because the same atomtype apparently can have different charges as well 
        void MakeResnameAtomNameToChargeMap();

        // make map from resname and atom name to atom type
        void MakeResnameAtomNameToTypeMap();

        // make map from residue name to atom type
        void MakeResidueToAtomTypeMap();

        // make the vector from indices to atomtypes
        void MakeIndicesToAtomType();

        Real getMassFromAtomTypeResname(const std::string& resname, const std::string& atomType);
        Real getChargeFromAtomNameResname(const std::string& resname, const std::string& atomName);
        std::string getAtomTypeFromAtomNameResname(const std::string& resname, const std::string& atomtype);

        Real getMassFromIndex(int index);
        Real getChargeFromIndex(int index);
        std::string getAtomTypeFromIndex(int index);

        enum TopIdx
        {
            atomtype = 1,
            resname = 3,
            atomName = 4,
            charge = 6,
            mass = 7
        };

        struct AtomType
        {
            using Real = CommonTypes::Real;

            std::string type_;
            std::string resname_;
            std::string atomName_;
            Real charge_;
            Real mass_;
        };

    private:
        std::string topName_;

        std::vector<std::string> resnames_;

        std::vector<AtomType> atomtypeIndices_;

        // Each line in [ atoms ] directive must have 8 entries
        int linenum_ = 8;

        std::vector<AtomType> atomtypes_;

        std::map<std::vector<std::string>, Real> ResNameAtomTypeToMassMap_;

        std::map<std::vector<std::string>, Real> ResNameAtomNameToChargeMap_;

        std::map<std::vector<std::string>, std::string> ResNameAtomNameToTypeMap_;

        std::map<std::string, int> MapResnameToNumber_;

        std::map<std::string, std::vector<AtomType>> MapResidueToAtomType_;
};