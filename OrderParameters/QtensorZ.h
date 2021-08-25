#include "Calculation.h"
#include "tools/Assert.h"
#include "xdr/MoleculeStructs.h"
#include "CalculationTools.h"

#include <vector>
#include <array>
#include <string>

// Class that calculates and average Qtensor as a function of z over time
// This can give P2 as a function of z(z can be any dimension)
// Or director as a function of z

class QtensorZ : public Calculation
{
    public:
        using Range = CommonTypes::Real2;

        QtensorZ(const CalculationInput& input);
        virtual ~QtensorZ(){};

        virtual void calculate() override;

    private:
        std::string residueName_;

        // The center of mass of each of the residues
        std::vector<Real3> COM_;

        // Get the Qtensor for each of the bins

        // The string that indicates the direction
        std::string direction_;
        int index_;

        // Map direction string to index
        std::map<std::string, int> MapNameToDirection = 
        {
            {"x", 0},
            {"y", 1},
            {"z", 2}
        };

        Range range_;
};