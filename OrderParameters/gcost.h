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
        void printHistogram2d(std::string name);
        void printrdfhist2d(std::string name);
        void printnumneighbors(std::string name);

        Real calcFactor(Real3& ui, Real3& uj);
        Real calcg1(Real3& ui, Real3& uj);
        Real calcg2(Real3& ui, Real3& uj);
        void registerCalcFunc(int i, fcn function);

        void initializeDistanceCOM();

    private:
        binptr bin_;
        binptr tbin_;
        int headindex1_=1;
        int tailindex1_=2;
        int headindex2_=1;
        int tailindex2_=2;

        int numatoms1_;
        int numatoms2_;
        int numresidues1_;
        int numresidues2_;

        // COM of whether or not a molecule is inside a PV
        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;
        std::vector<Real3> distanceCOM1_;
        std::vector<Real3> distanceCOM2_;

        // COMIndices for calculating the g(costheta)
        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;
        std::vector<int> distanceCOMIndices1_;
        std::vector<int> distanceCOMIndices2_;

        // volume for normalizing the r part 
        std::vector<Real> volume_;

        std::vector<Real3> uij1_;
        std::vector<Real3> uij2_;
        std::string residueName1_;
        std::string residueName2_;

        int numbins_;

        std::vector<int> InsideIndices1_;
        std::vector<int> InsideIndices2_;
        std::vector<int> OutsideIndices1_;
        std::vector<int> OutsideIndices2_;

        std::vector<Real> histogram_;
        std::vector<Real> histogramDotProduct_;
        std::vector<std::vector<Real>> histogramDotProduct2dPerIter_;
        std::vector<std::vector<Real>> histogramDotProduct2d_;
        std::vector<std::vector<Real>> histogramDotProduct2drdf_;

        std::map<int, fcn> MapIndexToFcn_;

        int index_=1;
        int numtbins_=20;

        Range trange_ = {{-1,1}};

        bool selfinteraction_=true;
};