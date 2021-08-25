#pragma once
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "xdr/GroFile.h"
#include "xdr/TopologyReader.h"
#include "AtomGroupParsingStrategy.h"
#include "xdr/MoleculeStructs.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

struct ResidueInput
{
    ParameterPack& pack_;
    GroFile& grofile_;
    TopologyReader& top_;
};

// ResidueGroup have the positions of all the atoms, the masses of all the masses and the atom indices of the atoms in that residue
class ResidueGroup
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using VectorReal3 = CommonTypes::VectorReal3;
        using stratptr = std::shared_ptr<AtomGroupParsingStrategy>;
        
        ResidueGroup(const ResidueInput& input);
        ~ResidueGroup(){};

        void update(const VectorReal3& total_atoms_);

        // getters
        std::string getName() const{return name_;}
        const std::vector<Molecule::residue> getResidues() const {return Residues_;}


    private:
        std::vector<int> AtomGroupGlobalIndices_;

        // Indices of the ResidueGroup
        std::vector<int> ResidueGroupIndices_;

        // number of atoms in the ResidueGroup
        int numResidueGroupatoms_;

        // name of the ResidueGroup
        std::string name_;

        // selection string of the ResidueGroup
        std::vector<std::string> selection_str_;

        // index string 
        std::vector<std::string> index_str_;

        // The residues
        std::vector<Molecule::residue> Residues_;

        // strategy for parsing
        stratptr strat_;

        //GroFile
        GroFile& grofile_;

        // The topology reader
        TopologyReader& top_;

        // sizes of each of the residues (number of atoms)
        std::vector<int> AtomsPerResidue_;

        // pointer to the atomgroup parsing strategy
        // mainly for resname and residue index
};