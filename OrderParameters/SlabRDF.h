#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "SimulationState.h"

#include <string>
#include <memory>

class SlabRDF : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;

        SlabRDF(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;
        void binUsingMinMax();

        void printSlabRDF(std::string name);

    private:
        std::string resname_;
        std::vector<std::vector<Real>> slabRDF_;
        std::vector<Real> slabN_;
        binptr zbin_;
        binptr rbin_;
        int numzbins_;
        int numrbins_;
        Real dz_;
        int headindex_=1;
        int tailindex_=2;
        Real above_=0;

        std::vector<Real> zBinLocation_;

        bool binUsingMinMax_=false;

        int direction_=2;
        std::vector<Real> volume_;

        bool within_z_=false;
};