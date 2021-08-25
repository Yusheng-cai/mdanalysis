#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "xdr/GroFile.h"
#include "xdr/TopologyReader.h"
#include "AtomGroupParsingStrategy.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

class ResidueGroup;

struct ResidueInput
{
    ParameterPack& pack_;
    GroFile& grofile_;
    TopologyReader& top_;
};

// AtomGroup is responsible for reading the input parameter pack and figuring out the correct indices for this particular AGroup
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

        // a residue would have information of all the atom positions as well as the masses & atomNumbers_
        struct residue
        {
            std::vector<Real3> atomPositions_;
            std::vector<Real>  atomMasses_;
            std::vector<int>   atomNumbers_;
        };

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
        std::vector<residue> Residues_;

        // strategy for parsing
        stratptr strat_;

        //GroFile
        GroFile& grofile_;

        // The topology reader
        TopologyReader& top_;

        // pointer to the atomgroup parsing strategy
        // mainly for resname and residue index
};