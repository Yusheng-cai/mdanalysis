#pragma once
#include "tools/Assert.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "SimulationState.h"
#include "xdr/TopologyReader.h"
#include "xdr/GroFile.h"
#include "AtomGroup.h"
#include "ResidueGroup.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <map>

struct CalculationInput
{
    ParameterPack& pack_;
    SimulationState& simstate_;
};

class Calculation
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        Calculation(const CalculationInput& input);
        virtual ~Calculation(){};

        virtual void update(){};
        virtual void calculate() = 0;

        void addAtomgroup(std::string name);
        void addResidueGroup(std::string name);

        const ResidueGroup& getResidueGroup(std::string name) const;

    protected:
        SimulationState& simstate_;

        std::map<std::string,int> MapAtomGroupNameToIndex_;
        std::map<std::string,int> MapResidueGroupNameToIndex_;

        std::vector<AtomGroup*> AtomGroups_; 
        std::vector<ResidueGroup*> ResidueGroups_;

        std::vector<std::string> ResidueNames_;
        std::vector<std::string> AtomGroupNames_;
};

namespace CalculationRegistry
{
    using Key = std::string;
    using Base= Calculation;

    using Factory = GenericFactory<Base, Key, const CalculationInput&>;

    template <typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const CalculationInput&>;
};