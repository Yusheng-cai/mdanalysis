#pragma once
#include "Calculation.h"
#include "tools/Assert.h"
#include "xdr/MoleculeStructs.h"
#include "CalculationTools.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "Qtensor.h"

#include <vector>
#include <iomanip>
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
        void binUsingMinMax();

        void printP2z(std::string name);
        void printPerIterP2z(std::ofstream& ofs);
        void printPerIterev(std::ofstream& ofs);
        void printPerIterNum(std::ofstream& ofs);
        void printPerIterQtensor(std::ofstream& ofs);

        void printPerIterP2zBeta(std::ofstream& ofs);
        void printPerItereveczBeta(std::ofstream& ofs);

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

        // record total number of residue per bin throughout simulation
        std::vector<Real> NumResPerBin_;

        std::vector<Matrix> BinnedMatrix_;

        // The head and tail index within the 5CB molecule
        int headIndex_;
        int tailIndex_;

        std::vector<Real> P2_;
        std::vector<Real> P2avg_;

        // Which atoms do you want to perform COM calculations over
        std::vector<int> COMIndex_;
        bool COMIndex_provided_=false;

        // ignore the P2 if the average N is less than this number
        Real ignoreP2LessThan_ = 0.0;

        int precision_ = 3;

        // nx^2, ny^2, nz^2 for the eigenvectors
        std::vector<Real3> eigvec_;

        // Per Iteration items
        std::vector<Matrix> BinnedMatrixIter_;
        std::vector<Real> NumResPerBinIter_;
        std::vector<Real> P2PerIter_;
        std::vector<Real3> evPerIter_;

        // zero Matrix
        Matrix zeroMatrix_ = {};

        // map from residue index to the corresponding bin index 
        std::vector<int> ResIndexToBinIndex_;

        int numbins_;
        bool usingMinMaxBins_=false;

        // we only bin the COM above this count
        Real above_=-100000;
};