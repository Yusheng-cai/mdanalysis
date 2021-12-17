#pragma once 

#include "Calculation.h"
#include "tools/Assert.h"
#include "SimulationState.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "Bin.h"

#include <string>
#include <memory>
#include <vector>
#include <array>

class DistributionGivenDimer : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;

        DistributionGivenDimer(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);

    private:
        Real rmax_;
        Real cosmax_;

        int headindex_;
        int tailindex_;

        std::string resName_;
        int numres_;

        std::vector<Real3> uij_;

        std::vector<Real> histogram_;

        binptr bin_;
        int numbins_;

        std::vector<Real3> COM_;

        // surface normal
        Real3 surfaceNormal_={{0,0,1}};

        // Angle with surface 
        std::vector<Real> AngleWithSurface;

        // Dimer per residue 
        std::vector<Real> DimerPerResidue_;
};