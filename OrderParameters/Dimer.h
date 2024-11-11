#pragma once

#include "Calculation.h"
#include "tools/CommonOperations.h"
#include "tools/CommonTypes.h"
#include "Bin.h"

#include <string>
#include <memory>
#include <vector>

class Dimer : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;

        Dimer(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        Real getNumDimer() {return dimer_vec_.size();}

        void printDimerIndices(std::ofstream& ofs);
        void printNonDimerIndices(std::ofstream& ofs);
        void printOrientation(std::string name);

    private:
        std::string resname_;
        Real alignment_cutoff_=-0.9;
        Real distance_cutoff_=0.35;
        Real distance_cutoff_B1B2_;

        int headindex_=1;
        int tailindex_=2;

        std::vector<Real3> uij_;
        int num_dimers_;
        std::vector<int> dimer_vec_;
        std::vector<int> non_dimer_vec_;
        std::vector<int> COMIndicesB1_={3,4,5,6,7,8};
        std::vector<int> COMIndicesB2_={9,10,11,12,13,14};

        std::vector<Real3> COMB1_;
        std::vector<Real3> COMB2_;

        std::vector<Real> orientation_monomer_;
        std::vector<Real> orientation_dimer_;
        int numtbins_=30;  
        Binptr bin_;
};