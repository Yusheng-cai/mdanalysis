#pragma once
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "parallel/OpenMP_buffer.h"
#include "xdr/GroFile.h"
#include "AtomGroupParsingStrategy.h"


#include <vector>
#include <memory>
#include <array>
#include <string>

class AtomGroup;

struct AtomGroupInput
{
    ParameterPack& pack_;
    GroFile& grofile_;
};

// AtomGroup is responsible for reading the input parameter pack and figuring out the correct indices for this particular AGroup
class AtomGroup
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using VectorReal3 = CommonTypes::VectorReal3;
        using stratptr  = AtomGroupParsingStrategy*; 
        
        AtomGroup(const AtomGroupInput& input);
        ~AtomGroup(){};

        void update(const VectorReal3& total_atoms_);

        // getters
        std::string getName() const{return name_;}
        const std::vector<int>& getAtomGroupIndices() const {return AtomGroupIndices_;} 
        const VectorReal3& getAtomPositions() const {return atom_positions_;}

        // accessors
        VectorReal3& accessAtomPositions() {return atom_positions_;} // perhaps this one should not be played around with

        // convert from AtomGroupIndices to global indices
        int AtomGroupIndices2GlobalIndices(int atomgroupIndices) const;


    private:
        // Indices of the AtomGroup
        std::vector<int> AtomGroupIndices_;

        // number of atoms in the AtomGroup
        int num_atoms_;

        // name of the AtomGroup
        std::string name_;

        // selection string of the AtomGroup
        std::vector<std::string> selection_str_;

        // index string 
        std::vector<std::string> index_str_;

        // positions of owned atoms
        VectorReal3 atom_positions_;

        // buffer for storing VectorReal3 when updating
        OpenMP::OpenMP_buffer<VectorReal3> atoms_positions_buffer_; 

        // Map Atom Group indices to global indices 
        std::map<int, int> AtomGroupIndicesToGlobalIndices_;

        // strategy for parsing
        stratptr strategy_;

        //GroFile
        GroFile& grofile_;
};