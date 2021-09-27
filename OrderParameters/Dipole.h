#include "Calculation.h"
#include "tools/Assert.h"
#include "Bin.h"
#include "LinAlgTools.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

class Dipole : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;

        Dipole(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
        void initializeBins();

    private:
        std::string residueName_;
        std::vector<Real3> dipoledirection_;
        std::vector<Real> histogram_;

        Real3 direction_;
        Binptr bin_;
};