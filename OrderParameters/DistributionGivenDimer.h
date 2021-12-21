#pragma once 

#include "Calculation.h"
#include "tools/Assert.h"
#include "SimulationState.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP.h"
#include "LinAlgTools.h"
#include "Bin.h"

#include <string>
#include <memory>
#include <vector>
#include <array>
/**
 * @brief Distribution of the cos(theta) given that an nCB molecule is r away from another molecule and forms an angle
 * cos(theta2) with a reference vector
 */

class DistributionGivenDimer : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;

        /**
         * @brief Construct a new Distribution Given Dimer object
         * 
         * @param input The calculation input pack 
         */
        DistributionGivenDimer(const CalculationInput& input);

        /**
         * @brief Perform calculation 
         * 
         */
        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
        void printNumDimerPerIter(std::ofstream& ofs);
        void printNumDimerPerResiduePerIter(std::ofstream& ofs);
        void printHistogramNotdimer(std::string name);
        void printDimerBetaFactor(std::ofstream& ofs);
        void printNotDimerBetaFactor(std::ofstream& ofs);

        Real getNumDimers() {return numDimersPerIter_;}

    private:
        Real rmax_;
        Real cosmax_;

        int headindex_=1;
        int tailindex_=2;

        std::string resName_;
        int numres_;
        // atoms per residue
        int numatoms_;

        std::vector<Real3> uij_;

        std::vector<Real> histogram_;
        std::vector<Real> histogramNotDimer_;

        binptr bin_;
        int numbins_;

        std::vector<Real3> COM_;
        std::vector<Real3> COMdistance_;

        // surface normal
        Real3 surfaceNormal_={{0,0,1}};

        // Angle with surface 
        std::vector<Real> AngleWithSurface_;

        // Dimer per residue 
        std::vector<Real> DimerPerResidue_;
        std::vector<Real> outsideDimerPerResidue_;

        // number of pairs of dimers formed per time step
        int numDimersPerIter_;

        std::vector<int> InsideIndices_;
        std::vector<int> OutsideIndices_;

        // indices of COM for calculating distances
        std::vector<int> DistanceCOM_;
};