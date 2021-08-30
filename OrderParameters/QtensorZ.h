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
        virtual void printOutput() override;

        std::vector<Real>& getP2() {return P2_;}

    private:
        std::string residueName_;

        // The center of mass of each of the residues
        std::vector<Real3> COM_;

        // Get the Qtensor for each of the bins

        // The string that indicates the direction
        std::string direction_ = "z";
        int index_;

        // name of the output
        std::string p2ZOutput_;
        std::ofstream p2zofs_;

        // Map direction string to index
        std::map<std::string, int> MapNameToDirection = 
        {
            {"x", 0},
            {"y", 1},
            {"z", 2}
        };

        Binptr bin_;

        // record total number of residue per bin throughout simulation
        std::vector<Real> NumResPerBin_;

        std::vector<Matrix> BinnedMatrix_;

        // The head and tail index within the 5CB molecule
        int headIndex_;
        int tailIndex_;

        // The indices of the atoms where COM is calculated wrt 
        std::vector<int> COMIndices_;

        std::vector<Real> P2_;

        // Which atoms do you want to perform COM calculations over
        std::vector<int> COMIndex_;
        bool COMIndex_provided_=false;

        // ignore the P2 if the average N is less than this number
        Real ignoreP2LessThan_ = 0.0;
};