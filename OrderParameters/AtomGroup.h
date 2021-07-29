#pragma once
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/Assert.h"

#include <vector>
#include <array>
#include <string>

// AtomGroup is responsible for reading the input parameter pack and figuring out the correct indices for this particular AGroup
class AtomGroup
{
    public:
        AtomGroup(const ParameterPack& pack);
        ~AtomGroup(){};

        void ParseSelectionString();

        // getters
        std::string getName() const{return name_;}
        const std::vector<int>& getAtomGroupIndices() const {return AtomGroupIndices_;} 


    private:
        // Indices of the AtomGroup
        std::vector<int> AtomGroupIndices_;

        // name of the AtomGroup
        std::string name_;

        // selection string of the AtomGroup
        std::vector<std::string> selection_str_;

        // index string 
        std::vector<std::string> index_str_;
};