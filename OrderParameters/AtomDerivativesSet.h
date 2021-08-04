#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <map>

struct AtomDerivatives
{
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;

    // This is the global index of the atom 
    int index;

    Real3 derivatives;
};

class AtomDerivativesSet
{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using VectorReal3 = CommonTypes::VectorReal3;

        AtomDerivativesSet(){};
        ~AtomDerivativesSet(){};

        // insert a derivative into AtomDerivatives
        void insert(int index, const Real3& derivative);
        void clear();

        const std::vector<AtomDerivatives>& getAtomDerivatives() const {return AtomDerivatives_;}
        std::vector<AtomDerivatives>& accessAtomDerivatives() {return AtomDerivatives_;};
        AtomDerivatives& accessAtomDerivativeByIndex(int index) {return AtomDerivatives_[index];}


    private:
        std::vector<AtomDerivatives> AtomDerivatives_;

        // an index map that keeps track of all the atom indices inserted
        std::map<int, bool> Index_Map_;
};