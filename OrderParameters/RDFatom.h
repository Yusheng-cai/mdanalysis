#include "Calculation.h"

#include <vector>
#include <array>
#include <string>

// calculate radial distribution function

// input can either be with residue or atom groups

// residuegroup1
// residuegroup2

// atomgroup1
// atomgroup2

class RDFatom : public Calculation
{
    public:
        RDFatom(const CalculationInput& input);

        virtual void calculate() override;

    private:
        std::string atomgroup1Name_;
        std::string atomgroup2Name_;
        std::string resgroup1Name_;
        std::string resgroup2Name_;
};