#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"
#include "Bin.h"
#include "SimulationState.h"

#include <string>
#include <memory>

class SlabResidueProperty : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Range2 = CommonTypes::Real2;

        SlabResidueProperty(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() {};
        virtual void finishCalculate() override;

        void binUsingMinMax(const std::vector<Real3>& positions);

        void printHistogram(std::string name);


    private:
        std::string resname_;
        Binptr bin_;
        Real dz_;
        int numbins_;
        Real above_;
        bool usingMinMax_=false;
        int direction_=2;
        std::vector<Real> ResidueLocationPerBin_;
        std::vector<Real> histogram_;
};