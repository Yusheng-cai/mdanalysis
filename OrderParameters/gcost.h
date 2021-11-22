#pragma once 
#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "parallel/OpenMP.h"
#include "LinAlgTools.h"

#include <vector>
#include <memory>
#include <array>
#include <string>

class gcost : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        gcost(const CalculationInput& input);

        virtual void calculate();
        virtual void finishCalculate() override;
        void printHistogram(std::string name);


    private:
        binptr bin_;
        int headindex_=1;
        int tailindex_=2;

        std::vector<Real3> COM_;
        std::vector<Real3> uij_;
        std::string residueName_;

        std::vector<std::vector<Real>> neighborDistance_;

        int numbins_;

        std::vector<int> InsideIndices_;
        OpenMP::OpenMP_buffer<std::vector<int>> InsideIndicesBuffer_;

        std::vector<std::vector<Real>> dotProduct_;

        std::vector<Real> histogram_;
        std::vector<Real> histogramDotProduct_;
        std::vector<int> histogramPerIter_;
        std::vector<Real> histogramDotProductPerIter_;
};