#include "Calculation.h"
#include "tools/Assert.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"

#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <numeric>
#include <memory>

class RDFcost : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>; 
        RDFcost(const CalculationInput& input);

        virtual void calculate() override;
        virtual void printOutput() override;

    private:
        std::string resName_;

        std::vector<int> COMIndices_;

        Binptr bin_;

        std::vector<Real> gcostheta_;

        std::vector<Real3> uij_;

        int headIndex_;
        int tailIndex_;
};