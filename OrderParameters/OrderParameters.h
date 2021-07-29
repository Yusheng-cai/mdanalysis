#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "Output_values.h"
#include "SimulationState.h"
#include "ProbeVolumeRegistry.h"

#include <map>
#include <string>
#include <functional>

struct OrderParametersInput
{
    const ParameterPack& pack_;
    SimulationState& simstate;
    const ProbeVolumeRegistry& pv_registry_;
};

class OrderParameters
{
    public:
        using OutputRegistry = std::map<std::string,OutputValue>; 
        using Real           = CommonTypes::Real;

        OrderParameters(const OrderParametersInput& input);
        virtual ~OrderParameters(){};

        virtual void calculate() = 0;
        virtual void update(){};

        // register the output value's function
        void registerOutput(std::string name, OutputValue::ValueFunction func);
        const OutputRegistry& getOutputRegistry() const {return OutputValues;} 

        // getters 
        Real getValue() const {return value;};

    protected:
        OutputRegistry OutputValues;

        Real value;
        
        SimulationState& simstate_;
};

namespace OrderParametersRegistry
{
    using Base = OrderParameters;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, const OrderParametersInput&>;
    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const OrderParametersInput&>;
}