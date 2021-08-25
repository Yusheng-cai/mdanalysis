#include "tools/Assert.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "OrderParameters/SimulationState.h"
#include "xdr/TopologyReader.h"
#include "xdr/GroFile.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>

struct CalculationInput
{
    ParameterPack& pack_;
    SimulationState& simstate_;
    TopologyReader& top;
    GroFile& gro;
};

class Calculation
{
    public:
        Calculation(const CalculationInput& input):top_(input.top), gro_(input.gro), simstate_(input.simstate_){};

        virtual void update(){};
        virtual void calculate() = 0;

    protected:
        TopologyReader& top_;
        GroFile& gro_;
        SimulationState& simstate_;
};

namespace CalculationRegistry
{
    using Key = std::string;
    using Base= Calculation;

    using Factory = GenericFactory<Base, Key, const CalculationInput&>;

    template <typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const CalculationInput&>;
};