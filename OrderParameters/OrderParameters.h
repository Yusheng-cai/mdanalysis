#pragma once
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "Output_values.h"
#include "SimulationState.h"
#include "AtomDerivativesSet.h"

#include <map>
#include <string>
#include <memory>
#include <functional>

struct OrderParametersInput
{
    const ParameterPack& pack_;
    SimulationState& simstate;
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
        const OutputValue& getOutput(std::string name_) const;
        void addAtomGroup(std::string name);

        // clear the derivatives stored in DerivativeOutputs
        void clearDerivativesOutputs();

        // clear a specific derivative for a specific atom group
        void clearDerivative(std::string& name);

        // get the index of the atom group
        int getAtomGroupIndex(std::string& name) const;

        // getters 
        Real getValue() const {return value;};
        std::string getName() const {return name_;}
        const std::vector<AtomDerivativesSet>& getDerivativeOutputs() const {return derivativeOutputs_;}
        std::vector<AtomDerivativesSet>& accessDerivativeOutputs() {return derivativeOutputs_;}

        const AtomDerivativesSet& getDerivatives(std::string& AtomGroupName) const; 
        AtomDerivativesSet& accessDerivatives(std::string& AtomGroupName);

        AtomGroup& accessAtomGroup(std::string& name); 
        const AtomGroup& getAtomGroup(std::string& name) const;

    protected:
        OutputRegistry OutputValues;

        Real value;
        
        SimulationState& simstate_;

        SimulationBox& simbox_;

        std::string name_;

        std::vector<AtomDerivativesSet> derivativeOutputs_;

        std::vector<const AtomGroup*> OPAtomGroups_;

        // a map from Atomgroup name to the index in the vector of atomderivativesSet
        std::map<std::string, int> MapAtomGroupNameToIndex_;
};

namespace OrderParametersRegistry
{
    using Base = OrderParameters;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, const OrderParametersInput&>;
    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const OrderParametersInput&>;
}