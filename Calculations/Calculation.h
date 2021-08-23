#include "tools/Assert.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>

struct CalculationInput
{
    ParameterPack& pack_;
};

class Calculation
{
    public:
        Calculation(const CalculationInput& input);

        virtual void update(){};
        virtual void calculate() = 0;

    protected:
};

namespace CalculationRegistry
{
    using Key = std::string;
    using Base= Calculation;

    using Factory = GenericFactory<Base, Key, const CalculationInput&>;

    template <typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const CalculationInput&>;
};