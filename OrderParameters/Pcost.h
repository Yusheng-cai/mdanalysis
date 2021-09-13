#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"

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

class Pcost : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;

        Pcost(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;
        virtual void printOutput() override;

    private:
        Binptr costBin_;

        // the indices for calculating the center of mass
        std::vector<int> COMIndices_;

        // the name of the residuegroup provided
        std::string residueGroupName_;

        // The center of mass
        std::vector<Real3> COM_;

        // histogram
        std::vector<Real> histogram_;

        int headIndex_;
        int tailIndex_;

        int directionIndex_;

        std::vector<Real3> uij_;

        // the vector in which all the molecules are to calculate their cos theta against
        // assumes z direction if not specified
        std::array<Real,3> arr_ = {{0,0,1}};

        // output stream
        std::ofstream ofs_;
        std::string outputName_;
        int precision_ = 3;

        // number of residues per bin
        std::vector<Real> numResiduePerBin_;

        // name of probe volume
        std::string ProbeVolumeName_;
};