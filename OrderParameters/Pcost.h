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

        void printHistogram(std::string name);
        void printAtomIndices(std::string name);

        void printcosthetaHistogramPerIter(std::ofstream& ofs);

    protected:
        Binptr costBin_;

        // the name of the residuegroup provided
        std::string residueGroupName_;

        // histogram
        std::vector<Real> histogram_;

        // histogram per iteration
        std::vector<Real> histogramPerIter_;

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

        // name of probe volume
        std::string ProbeVolumeName_;

        // name of the atom Indices output
        std::string AtomIndicesName_;

        // keeps track of AtomIndicesInPV per frame
        std::vector<std::vector<int>> AtomIndicesInPV_;

        std::vector<int> InsideIndices_;
};