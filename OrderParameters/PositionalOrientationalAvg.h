#pragma once 

#include "Calculation.h"
#include "LinAlgTools.h"
#include "ResidueGroup.h"
#include "SimulationState.h"
#include "Bin.h"

#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <map>

/*
    Class that calculates the positional orientational average

    Mathematical this is represented by 
        <O(x) | CosTheta, R> = 1/N_ij \int O(x) \delta(Cos(\theta)_i -  CosTheta) \delta (r_ij - R) 
*/
class PositionalOrientationalAvg : public  Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using avg_func = std::function<Real(const Molecule::residue&, const Molecule::residue&)>; 

        PositionalOrientationalAvg(const CalculationInput& input);

        void RegisterOutputs();
        void RegisterCalculationFunctions(std::string name, avg_func func);

        virtual void calculate();
        virtual void update(){};
        virtual void finishCalculate();

        // calculation functions 
        Real PairUsr(const Molecule::residue& res1, const Molecule::residue& res2);

        // printing outputs
        void printAverage(std::string name);

    private:
        // The indices inside a residue to calculate Usr  
        std::vector<int> UsrIndices_;

        // The name of the residue 
        std::string ResidueName_;

        // Histogram
        std::vector<std::vector<Real>> PosOriHist_;
        std::vector<std::vector<Real>> PosOriAvg_;

        // Bins
        Binptr RBin_;
        Binptr CosThetaBin_;
        int NumRBins_;
        int NumCosThetaBins_;

        // director of each of the molecule
        std::vector<Real3> uij_;

        // Head and Tail index
        int HeadIndex_=1;
        int TailIndex_=2;

        // The kind of averaging that user wants to perform
        std::string Avg_type_;
        std::map<std::string, avg_func> MapNameToFunc_;

        // factors for Usr
        Real Usr_factor_ = 138.935458; 
        Real Usr_cuttoff_= 1.2;
        Real Usr_epsilon_= 0.384;
        Real Usr_alpha_;
};