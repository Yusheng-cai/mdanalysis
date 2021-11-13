#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"
#include "BetaFactorWriter.h"

#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>

// calculate distribution of cosine theta of the angles in a probe volume

class Cost : public Calculation
{
    public:
        using BFptr = std::unique_ptr<BetaFactorWriter>;
        using Binptr = std::unique_ptr<Bin>;


        Cost(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;
        virtual void printOutputOnStep() override;
        void printhistogram(std::string name);
        void printavgCosthetaPerIter(std::ofstream& ofs);
        void printNtilde(std::ofstream& ofs);

        Real getCostheta2() {return avgCostheta_;}

    private:
        // bin pointer
        Binptr Bin_;
        int numBins_ = 50;

        // the indices for calculating the center of mass
        std::vector<int> COMIndices_;

        // the name of the residuegroup provided
        std::string residueGroupName_;

        // The center of mass
        std::vector<Real3> COM_;

        int headIndex_;
        int tailIndex_;

        int directionIndex_;

        std::vector<Real3> uij_;

        // the vector in which all the molecules are to calculate their cos theta against
        // assumes z direction if not specified
        std::array<Real,3> arr_ = {{0,0,1}};

        // output Atom Indices
        std::ofstream ofs_;

        // name of probe volume
        std::string ProbeVolumeName_;

        // name of the atom Indices output
        std::string BetaFactorsName_;

        // keeps track of AtomIndicesInPV per frame
        std::vector<Real> BetaFactors_;

        int precision_ = 3;

        // histogram 
        std::vector<Real> histogram_;

        BFptr bf_;

        Real avgCostheta_;
};