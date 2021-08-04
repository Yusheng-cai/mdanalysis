#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP_buffer.h"

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

        std::size_t size() const {return AtomDerivatives_.size();}

        const std::vector<AtomDerivatives>& getAtomDerivatives() const {return AtomDerivatives_;}
        const AtomDerivatives& getAtomDerivativeByIndex(int index) const {return AtomDerivatives_[index];} 
        std::vector<AtomDerivatives>& accessAtomDerivatives() {return AtomDerivatives_;};
        AtomDerivatives& accessAtomDerivativeByIndex(int index) {return AtomDerivatives_[index];}

        void insertOMP(int index, const Real3& derivative);
        void CombineAndClearOMPBuffer();
        void setMasterObject();


    private:
        std::vector<AtomDerivatives> AtomDerivatives_;

        // This is used solely for OMP regions
        OpenMP::OpenMP_buffer<std::vector<AtomDerivatives>> AtomDerivativesBuffer_;

        // an index map that keeps track of all the atom indices inserted
        std::map<int, bool> Index_Map_;
};