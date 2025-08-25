#pragma once 
#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "Bin.h"
#include "parallel/OpenMP.h"
#include "LinAlgTools.h"
#include "SimulationState.h"
#include "tools/Constants.h"

#include <vector>
#include <memory>
#include <array>
#include <string>
#include <functional>
#include <map>

class DipoleDistribution : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        using fcn    = std::function<Real(Real3&, Real3&)>;
        using range  = std::array<Real,2>;

        DipoleDistribution(const CalculationInput& input);

        virtual void calculate();
        virtual void finishCalculate() override;

        void printHistogram2d(std::string name);

    private:
        binptr bin_;
        binptr tbin_;
        int headindex_=1;
        int tailindex_=2;

        int numatoms_;
        int numresidues_;

        // COM of whether or not a molecule is inside a PV
        std::vector<Real3> COM_;
        std::vector<int> COMIndices_;
        std::vector<Real3> uij_;
        std::string residueName_;
        range trange_={0,1};

        int numbins_;

        // I want the output to be (numrbins, numtbins)
        std::vector<std::vector<Real>> histogram2d_;

        // calculation functions
        int numtbins_=50;
};