#pragma once 

#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include "Qtensor.h"
#include "Bin.h"
#include "parallel/OpenMP.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

class grcost : public Calculation 
{
    public:
        using binptr = std::unique_ptr<Bin>;
        using Matrix = CommonTypes::Matrix;
        grcost(const CalculationInput& input);

        virtual void calculate() override;

    private:
        std::vector<Real3> COM_;
        std::vector<Real3> uij_;

        binptr tbin_;
        binptr rbin_;

        std::string residueName_;

        int headindex_=1;
        int tailindex_=2;

        Matrix Qtensor_;

        std::vector<int> InsideIndices_;

        Real3 director_;

        std::vector<std::vector<Real>> rdotdirector_;
        std::vector<std::vector<Real>> dotproduct_;
        std::vector<std::vector<int>> histogramPerIter_;
        std::vector<std::vector<Real>> histogramDotProductPerIter_;
        std::vector<std::vector<Real>> neighborDistance_;
};