#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "SimulationState.h"
#include "LinAlgTools.h"

#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>

/*
    This class calculate the distribution function P(z,cos(beta))=<d(z-z_i) * d(cos(\beta) - u_i * z_i)>

    definition:    
    z is the height of the LC molecule in the slab
    cos(beta) is the angle that the LC forms with respect to some user defined direction.
*/
class SlabOrientation : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Range2 = std::array<Real,2>;

        SlabOrientation(const CalculationInput& input);

        /*
            @Function that registers outputs
        */
        void RegisterOutputs();

        /*
            Function that reads the inputs 
        */
        void ReadInputs();

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
        void printNumResidue(std::string name);
        void printHistogramSquared(std::string name);

        void binUsingMinMax();

    private:
        Binptr costBin_;
        Binptr costsquaredBin_;
        Binptr zBin_;

        std::string direction_ = "z";

        std::map<std::string, int> MapdirectionToIndex_ = {
            {"x", 0},
            {"y", 1}, 
            {"z", 2}
        };

        // the name of the residuegroup provided
        std::string residueGroupName_;

        // This is 2d histogram for P(cos(theta), z)
        std::vector<std::vector<Real>> histogram2d_;
        std::vector<std::vector<Real>> BWcostZ_;

        // This is 2d histogram for P(cos2(theta),z)
        std::vector<std::vector<Real>> histogram2d_squared_;

        int headIndex_;
        int tailIndex_;

        int directionIndex_;

        std::vector<Real3> uij_;

        // the vector in which all the molecules are to calculate their cos theta against
        // assumes z direction if not specified
        std::array<Real,3> arr_ = {{0,0,1}};

        // output stream
        int precision_ = 5;

        // number of residues per bin
        std::vector<Real> numResiduePerBin_;
        std::vector<Real> numResiduePerBinSquared_;
        std::vector<Real> ResidueLocationPerBin_;

        // we are binning z directions using min/max of the COM
        int numzbins_;
        int numtbins_=30;
        Real above_;
        bool usingMinMax_=false;

        // anchoring strenght 
        std::vector<Real> anchoring_strength_;
};