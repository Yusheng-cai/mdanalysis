#pragma once
#include "Calculation.h"
#include "tools/Assert.h"
#include "xdr/MoleculeStructs.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "SimulationState.h"
#include "LinAlgTools.h"

#include <vector>
#include <iomanip>
#include <array>
#include <string>
#include <memory>

// Class that calculates and average Qtensor as a function of z over time
// This can give P2 as a function of z(z can be any dimension)
// Or director as a function of z
class SlabQtensor : public Calculation
{
    public:
        using Range = CommonTypes::Real2;
        using Binptr= std::unique_ptr<Bin>;
        using Matrix= CommonTypes::Matrix;

        SlabQtensor(const CalculationInput& input);
        virtual ~SlabQtensor(){};

        virtual void calculate() override;
        virtual void finishCalculate() override;
        void binUsingMinMax();

        void printP2z(std::string name);
        void printP2zbeta(std::string name);
        void printevBeta(std::string name);

        std::vector<Real>& getP2() {return P2_;}

    private:
        std::string residueName_;

        // The center of mass of each of the residues
        std::vector<Real3> COM_;

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

        // record total number of residue per bin throughout simulation
        std::vector<Real> NumResPerBin_;

        // Qtensor as a function of z 
        std::vector<Matrix> QtensorZ_;

        // The head and tail index in the 5CB molecule
        int headIndex_=1;
        int tailIndex_=2;

        std::vector<Real> P2_;
        std::vector<Real> P2avg_;

        int precision_ = 3;

        // nx^2, ny^2, nz^2 for the eigenvectors
        std::vector<Real3> eigvec_;

        // map from residue index to the corresponding bin index 
        std::vector<int> ResIndexToBinIndex_;

        // number of bins if we are binning it from min to max 
        int numbins_;
        bool usingMinMaxBins_=false;
        std::vector<Real> BinLocation_;

        // we only bin the COM above this count
        Real above_=-100000;
};