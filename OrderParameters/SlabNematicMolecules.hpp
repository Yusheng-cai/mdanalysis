#include "Calculation.h"
#include "Bin.h"
#include "tools/CommonTypes.h"
#include "IndicatorFunction1d.h"

#include <memory>

class SlabNematicMolecules : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Range2 = CommonTypes::Real2;

        SlabNematicMolecules(const CalculationInput& input);
        void ReadInputs();

        virtual void calculate();

        void binUsingMinMax();

        virtual void finishCalculate();

        void printAverageVectorStilde(std::string inputfname);
    private:
        IndicatorFunction1d func_left_;
        IndicatorFunction1d func_right_;

        int headIndex_=1;
        int tailIndex_=2;

        int directionIndex_;

        std::vector<Real3> uij_;

        // the vector in which all the molecules are to calculate their cos theta against
        // assumes z direction if not specified
        std::array<Real,3> arr_ = {{0,0,1}};

        // output stream
        int precision_ = 5;

        // number of residues per bin
        std::vector<Real> numResiduePerBin_;
        std::vector<Real> numResiduePerBinSquared_;
        std::vector<Real> ResidueLocationPerBin_;

        // we are binning z directions using min/max of the COM
        int numzbins_;
        int numtbins_=30;
        Real above_;
        bool usingMinMax_=false;

        Binptr zBin_;
        std::string residueGroupName_;
        std::string direction_ = "z";

        std::map<std::string, int> MapdirectionToIndex_ = {
            {"x", 0},
            {"y", 1}, 
            {"z", 2}
        };

        Real left_max_=-0.7, right_max_=0.7;
        Real left_sigma_=0.01, right_sigma_=0.01;
        Real left_alphaC_ = 0.02, right_alphaC_=0.02;

        std::vector<Real> vector_stilde_;
        std::vector<Real> vector_stilde_iter_;
};