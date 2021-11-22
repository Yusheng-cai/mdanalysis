#pragma once 

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Calculation.h"
#include "tools/InputParser.h"
#include "SimulationState.h"

#include <vector>
#include <string>

class findNearStructure : public Calculation
{
    public:
        findNearStructure(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override {};

        void calculateAtom();
        void calculateRes();

        void printAtomIndicesPerIter(std::ofstream& ofs);

    private:
        std::string agName_;
        std::string pvName_;

        // residue name
        std::string resName_;

        // The COM of the residues 
        std::vector<Real3> COM_;

        // The indices of the atoms near structure 
        std::vector<int> AtomIndices_;

        std::string mode_ = "residue";
};