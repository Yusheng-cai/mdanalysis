#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Calculation.h"
#include "tools/InputParser.h"
#include "CalculationTools.h"

#include <vector>
#include <string>

class findNearStructure : public Calculation
{
    public:
        findNearStructure(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override {};

        void printAtomIndicesPerIter(std::ofstream& ofs);

    private:
        std::string resname_;
        std::string pvName_;

        // The COM of the residues 
        std::vector<Real3> COM_;

        // The indices of the atoms near structure 
        std::vector<int> AtomIndices_;
};