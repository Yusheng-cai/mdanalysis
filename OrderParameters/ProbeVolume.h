#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/GenericFactory.h"
#include "SimulationState.h"
#include "SimulationBox.h"

#include <vector>
#include <array>
#include <string>

struct ProbeVolumeOutput
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using VectorReal3 = CommonTypes::VectorReal3;

    Real hx_;
    Real htilde_x_;
    Real3 dhtilde_dx_;
    VectorReal3 dyn_derivatives_;
};

struct ProbeVolumeInput
{
    ParameterPack& ParamPack;
    SimulationState& simstate;
};

class ProbeVolume
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        // input for ProbeVolume 
        ProbeVolume(ProbeVolumeInput& input):simstate_(input.simstate), simbox_(simstate_.getSimulationBox()){};
        virtual ~ProbeVolume(){};

        // Update the ProbeVolume as needed -> usually used for Dynamic Probe Volumes
        virtual void update(){};
        virtual ProbeVolumeOutput calculate(const Real3& x) = 0;

        bool isDynamic(){return isDynamic_;}

    private:
        // default Dynamic to be false
        bool isDynamic_=false;
        SimulationState& simstate_;
        SimulationBox& simbox_;
};

namespace ProbeVolumeRegistry
{
    using Base    = ProbeVolume;
    using Key     = std::string;

    using Factory = GenericFactory<Base,Key, ProbeVolumeInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, ProbeVolumeInput&>;
}