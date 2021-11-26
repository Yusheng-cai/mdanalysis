#pragma once 
#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "Bin.h"
#include "parallel/OpenMP.h"
#include "LinAlgTools.h"
#include "SimulationState.h"

#include <vector>
#include <memory>
#include <array>
#include <string>
#include <functional>
#include <map>

class gcost : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        using fcn    = std::function<Real(Real3&, Real3&)>;
        using range  = std::array<Real,2>;
        gcost(const CalculationInput& input);

        virtual void calculate();
        virtual void finishCalculate() override;
        void printHistogram(std::string name);

        Real calcFactor(Real3& ui, Real3& uj);
        Real calcg1(Real3& ui, Real3& uj);
        Real calcg2(Real3& ui, Real3& uj);
        void registerCalcFunc(int i, fcn function);

    private:
        binptr bin_;
        binptr tbin_;
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

        std::vector<std::vector<Real>> histogramDotProduct2d_;
        OpenMP::OpenMP_buffer<std::vector<std::vector<Real>>> histogramDotProduct2dbuffer_;

        OpenMP::OpenMP_buffer<std::vector<Real>> histogramDotProductPerIterbuffer_;
        OpenMP::OpenMP_buffer<std::vector<int>> histogramPerIterbuffer_;

        std::map<int, fcn> MapIndexToFcn_;

        int index_=1;
        int numtbins_=20;

        Range trange_ = {{-1,1}};
};