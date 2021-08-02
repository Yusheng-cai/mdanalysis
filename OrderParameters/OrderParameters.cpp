#include "OrderParameters.h"

OrderParameters::OrderParameters(const OrderParametersInput& input)
:simstate_(input.simstate), simbox_(input.simstate.getSimulationBox())
{}

void OrderParameters::registerOutput(std::string name, OutputValue::ValueFunction func)
{
    OutputValue output(name, func);

    OutputValues.insert(std::make_pair(name, output)); 
}

const OutputValue& OrderParameters::getOutput(std::string name) const
{
    auto it = OutputValues.find(name);

    ASSERT((it != OutputValues.end()), "The name " << name << " is not available.");

    return it -> second;
}