#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "Qtensor.h"
#include "SimulationState.h"

#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>

class Pcostz : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Range2 = std::array<Real,2>;

        Pcostz(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
        void printHistogramPerIter(std::ofstream& ofs);

        void binUsingMinMax();

    private:
        Binptr costBin_;
        Binptr zBin_;

        std::string direction_ = "z";

        std::map<std::string, int> MapdirectionToIndex_ = {
            {"x", 0},
            {"y", 1}, 
            {"z", 2}
        };

        // the name of the residuegroup provided
        std::string residueGroupName_;

        std::vector<std::vector<Real>> histogram2d_;
        std::vector<std::vector<Real>> histogramIter_;

        int headIndex_;
        int tailIndex_;

        int directionIndex_;

        std::vector<Real3> uij_;

        // the vector in which all the molecules are to calculate their cos theta against
        // assumes z direction if not specified
        std::array<Real,3> arr_ = {{0,0,1}};

        // output stream
        int precision_ = 3;

        // number of residues per bin
        std::vector<Real> numResiduePerBin_;

        // number of residues per bin per iteration
        std::vector<Real> numResiduePerBinIter_;

        int ignoreBelow_ = 0;

        // we are binning z directions using min/max of the COM
        int Znumbins_;
        Real above_;
        bool usingMinMax_=false;
};