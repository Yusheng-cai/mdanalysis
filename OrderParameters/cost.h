#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "BetaFactorWriter.h"
#include "LinAlgTools.h"
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

// calculate distribution of cosine theta of the angles in a probe volume

class Cost : public Calculation
{
    public:
        using BFptr = std::unique_ptr<BetaFactorWriter>;
        using Binptr = std::unique_ptr<Bin>;

        Cost(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printhistogram(std::string name);
        void printavgCosthetaPerIter(std::ofstream& ofs);
        void printNtilde(std::ofstream& ofs);

        Real getavgCostheta2() const {return avgCostheta_;}
        Real getNumCOM() const {return numCOM_;}
        Real getCostheta2() const {return Costheta_;}

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

        // name of probe volume
        std::string ProbeVolumeName_;

        int precision_ = 3;

        // histogram 
        std::vector<Real> histogram_;
        Real avgCostheta_;
        Real Costheta_;

        int numCOM_;
};