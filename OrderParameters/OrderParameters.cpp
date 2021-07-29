#include "OrderParameters.h"

OrderParameters::OrderParameters(const OrderParametersInput& input)
:simstate_(input.simstate)
{}

void OrderParameters::registerOutput(std::string name, OutputValue::ValueFunction func)
{
    OutputValue output(name, func);

    OutputValues.insert(std::make_pair(name, output)); 
}