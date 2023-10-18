#include "Calculation.h"
#include "IndicatorFunction1d.h"
#include "tools/CommonTypes.h"

#include <string>

class NematicMolecules: public Calculation
{
    public:
        NematicMolecules(const CalculationInput& input);
        virtual void calculate() override;
        virtual void update() {};
        virtual void finishCalculate() {};

        Real get_stilde() {return stilde_;}
    private:
        IndicatorFunction1d func_left_;
        IndicatorFunction1d func_right_;
        Real3 arr_;

        std::string ResidueGroupName_;
        int HeadIndex_, TailIndex_;
        bool useDirector_=false;

        std::vector<Real> stilde_i_;

        std::vector<Real3> uij_;

        Real stilde_;

        std::vector<int> AtomIndices_;

        Real left_max_=-0.7, right_max_=0.7;
        Real left_sigma_=0.01, right_sigma_=0.01;
        Real left_alphaC_ = 0.02, right_alphaC_=0.02;
};