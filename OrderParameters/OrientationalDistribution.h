#pragma once

#include "Calculation.h"
#include "Bin.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "SimulationState.h"
#include "ResidueGroup.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

/*
    Class that finds the distribution of LC molecules within a certain region of interest (probe volume)
*/
class OrientationalDistribution : public Calculation
{
    public:
        using Binptr  = std::unique_ptr<Bin>;
        using Range   = std::array<Real,2>;
        using Matrix  = CommonTypes::Matrix;

        OrientationalDistribution(const CalculationInput& input);

        void registerOutputs();
        void registerOutputfile();

        virtual ~OrientationalDistribution(){};
        virtual void calculate();
        virtual void update();
        virtual void finishCalculate() override;

        // printing functions
        void PrintDistribution(std::string name);
        void PrintCosthetasquared_betafactors(std::ofstream& ofs);
        void PrintDistributionXY(std::string name);

        // per iterations printing functions
        void printResidueAngles(std::ofstream& ofs);

        // printing function for output files 
        Real getAvgCostheta() {return AvgCostheta_;}
        Real getAvgCosthetasquared() {return AvgCosthetasquared_;}

    private:
        std::vector<Real> PCosTheta_;
        std::vector<Real> PCosThetaSquared_;

        // head and tail index used for identifying orientation of an LC molecule 
        int HeadIndex_=1;
        int TailIndex_=2;

        // name of the residue group
        std::string ResidueGroupName_;

        // number of bins
        int NumBins_;

        Real AvgCostheta_=0.0;
        Real AvgCosthetasquared_=0.0;

        // bins 
        Binptr CosThetaBin_;
        Range  CosThetaRange_ = {{-1.0,1.0}};
        Binptr CosThetaSquaredBin_;
        Range CosThetaSquaredRange_ = {{0.0,1.0}};

        // distribution xy 
        std::vector<std::vector<Real>> PcostDistributionXY_;
        std::vector<std::vector<Real>> HisotgramXY_;
        int numxbin_=50;
        int numybin_=50;
        Binptr xbin_;
        Binptr ybin_;

        // director of each of the molecules
        std::vector<Real3> uij_;

        // array that we want to calculate cos theta on
        Real3 arr_={{0,0,1}};
        bool useDirector_=true;

        // these are the costhetasquared of each of the atom for the residues in question
        std::vector<Real> costhetasquared_betafactor_;

        std::vector<int> AtomIndices_;

        std::map<int, Real> MapIndexToAngle_;
};