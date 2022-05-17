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

struct AtomGroupInput
{
    ParameterPack& pack_;
    GroFile& grofile_;
};

namespace OP
{
    struct Atom
    {
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        Real3 position;
        int index;
    };
};

// AtomGroup is responsible for reading the input parameter pack and figuring out the correct indices for this particular AGroup
class AtomGroup
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using VectorReal3 = CommonTypes::VectorReal3;
        using stratptr  = std::shared_ptr<AtomGroupParsingStrategy>; 
        
        AtomGroup(const AtomGroupInput& input);
        ~AtomGroup(){};

        void update(const VectorReal3& total_atoms_);

        // getters
        std::string getName() const{return name_;}
        const std::vector<int>& getAtomGroupIndices() const {return AtomGroupGlobalIndices_;} 
        int getNumAtoms() const {return atoms_.size();}

        // convert from AtomGroupIndices to global indices
        int AtomGroupIndices2GlobalIndices(int atomgroupIndices) const;
        int GlobalIndices2AtomGroupIndices(int globalIndices) const;

        const OP::Atom getAtomByIndex(int index) const {return atoms_[index];} 
        OP::Atom accessAtomByIndex(int index) {return atoms_[index];}
        const std::vector<OP::Atom>& getAtoms() const {return atoms_;}
        std::vector<OP::Atom>& accessAtoms() {return atoms_;}

        

    private:
        // Indices of the AtomGroup
        std::vector<int> AtomGroupGlobalIndices_;

        // number of atoms in the AtomGroup
        int numAtomGroupatoms_;

        // name of the AtomGroup
        std::string name_;

        // selection string of the AtomGroup
        std::vector<std::string> selection_str_;

        // index string 
        std::vector<std::string> index_str_;

        // The atoms
        std::vector<OP::Atom> atoms_;

        OpenMP::OpenMP_buffer<std::vector<OP::Atom>> atoms_buffer_;

        // Map Atom Group indices to global indices 
        std::map<int, int> AtomGroupIndicesToGlobalIndices_;
        std::map<int, int> GlobalIndicesToAtomGroupIndices_;

        // strategy for parsing
        stratptr strategy_;

        //GroFile
        GroFile& grofile_;
};