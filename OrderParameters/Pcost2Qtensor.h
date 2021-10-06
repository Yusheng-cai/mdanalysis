#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"
#include "Pcost.h"
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

class Pcost2Qtensor : public Calculation
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real3  = CommonTypes::Real3;
        using Real   = CommonTypes::Real;
        using bptr   = std::unique_ptr<BetaFactorWriter>; 

        Pcost2Qtensor(const CalculationInput& input);

        void printcost20PerResIter(std::ofstream& ofs);
        void printcost21PerResIter(std::ofstream& ofs);
        void printcost22PerResIter(std::ofstream& ofs);
        virtual void calculate() override;
        virtual void finishCalculate() override {};
        virtual void printOutputOnStep();
    
    private:
        Matrix Qtensor_;

        // The director of each of the molecules
        std::vector<Real3> uij_;
        std::vector<Real> norm_;
        std::vector<std::vector<Real>> cost2Peratom_;
        std::vector<std::vector<Real>> cost2PerRes_;

        // The eigenvector directions 
        Matrix eigenvector_;
        Real3 v0_;
        Real3 v1_;
        Real3 v2_;

        int head_index_;
        int tail_index_;

        std::string residue_name_;

        int atomSize_;
        int size_;

        bptr bf_;
};