#pragma once 

#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "CalculationTools.h"
#include "LinAlgTools.h"
#include "Qtensor.h"
#include "parallel/OpenMP.h"

#include <iostream>
#include <vector>
#include <array>
#include <string>

class LocalOrder : public Calculation
{
    public:
        using Matrix = CommonTypes::Matrix;

        LocalOrder(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override {};

        void printLocalOrderBetaFactor(std::ofstream& ofs);
        void printLocalBiaxialityBetaFactor(std::ofstream& ofs);
        void printLocaldirector(std::ofstream& ofs);

    private:
        Real radius_;

        std::string residueName_;

        // head index 
        int headIndex_;
        int tailIndex_;

        std::vector<Real3> uij_;

        // array is always default to be the z axis 
        Real3 arr_={{0,0,1}};

        std::vector<std::vector<int>> neighborIndices_;

        std::vector<Real> localP2_;
        std::vector<Real> localBiaxiality_;
        std::vector<Real> localdirector_;

        OpenMP::OpenMP_buffer<std::vector<std::vector<int>>> NeighborIndicesBuffer_;
};