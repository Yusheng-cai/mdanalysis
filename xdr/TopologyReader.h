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

        // make map from atomName to type
        void MakeAtomNameToTypeMap();

        // make map from atomName to mass
        void MakeAtomNameToMassMap();

        // make map from atomName to charge
        void MakeAtomNameToChargeMap();

        Real getMassFromAtomName(const std::string& atomName);
        Real getChargeFromAtomName(const std::string& atomName);

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

        // Each line in [ atoms ] directive must have 8 entries
        int linenum_ = 8;

        std::vector<AtomType> atomtypes_;

        std::map<std::string, std::string> AtomNameToTypeMap_;

        std::map<std::string, Real> AtomNameToMassMap_;

        std::map<std::string, Real> AtomNameToChargeMap_;
};