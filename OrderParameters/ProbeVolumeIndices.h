#pragma once 

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Calculation.h"
#include "tools/InputParser.h"
#include "SimulationState.h"

#include <vector>
#include <string>

class ProbeVolumeIndices : public Calculation
{
    public:
        ProbeVolumeIndices(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override {};

        void calculateAtom();
        void calculateRes();
        // This function assumes that a residue is inside the PV if any of its atoms is inside the PV 
        void calculateAtomRes();

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