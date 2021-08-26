#pragma once
#include "Calculation.h"
#include "tools/Assert.h"
#include "xdr/MoleculeStructs.h"
#include "CalculationTools.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "Qtensor.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

// Class that calculates and average Qtensor as a function of z over time
// This can give P2 as a function of z(z can be any dimension)
// Or director as a function of z

class QtensorZ : public Calculation
{
    public:
        using Range = CommonTypes::Real2;
        using Binptr= std::unique_ptr<Bin>;
        using Matrix= CommonTypes::Matrix;

        QtensorZ(const CalculationInput& input);
        virtual ~QtensorZ(){};

        virtual void calculate() override;
        virtual void finishCalculate() override;

        std::vector<Real>& getP2() {return P2_;}

    private:
        std::string residueName_;

        // The center of mass of each of the residues
        std::vector<Real3> COM_;

        // Get the Qtensor for each of the bins

        // The string that indicates the direction
        std::string direction_ = "z";
        int index_;

        // Map direction string to index
        std::map<std::string, int> MapNameToDirection = 
        {
            {"x", 0},
            {"y", 1},
            {"z", 2}
        };

        Binptr bin_;

        std::vector<Matrix> BinnedMatrix_;

        // The head and tail index within the 5CB molecule
        int headIndex_;
        int tailIndex_;

        // The indices of the atoms where COM is calculated wrt 
        std::vector<int> COMIndices_;

        std::vector<Real> P2_;
};