#pragma once

#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "xdr/MoleculeStructs.h"
#include "Bin.h"
#include "SimulationState.h"
#include "LinAlgTools.h"

#include <vector>
#include <array>
#include <memory>
#include <string>

class SRE_dimer : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>; 
        using Int2   = std::array<int,2>;

        SRE_dimer(const CalculationInput& input);
        ~SRE_dimer(){};

        virtual void calculate() override;
        void calculateSRE(int i, int j, Real& attr, Real& repul, Real& total);

        virtual void finishCalculate() override;

        void printHistograms(std::string name);

    private:
        Binptr bin_;
        std::string residueName_; 

        int head_index_=1;
        int tail_index_=2;

        Real Rlim_;
        Real Rlimsq_;

        // number of bins 
        int binnum_;

        // histogram for attractive , repulsive , total part of SRE
        std::vector<Real> histogram_attr_;
        std::vector<Real> histogram_repu_;
        std::vector<Real> histogram_tota_;

        // count for each of the bins
        std::vector<Real> count_;
        std::vector<Real3> director_;

        // constants for SRE
        Real factor_=138.9354;
        Real sigma_=0.384;
        Real beta_;
};