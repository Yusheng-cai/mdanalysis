#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/GenericFactory.h"
#include "SimulationState.h"
#include "SimulationBox.h"

#include <vector>
#include <array>
#include <iostream>
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
    const ParameterPack& ParamPack;
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
        virtual ProbeVolumeOutput calculate(const Real3& x) const = 0;
        virtual void setGeometry(){};

        bool isDynamic(){return isDynamic_;}

        // Obtain the simulation state reference
        const SimulationState& getSimulationState() const{return simstate_;}
        const SimulationBox& getSimulationBox() const{return simbox_;}

        // getters
        Real getSigma() const {return sigma_;}
        Real getAlphaC() const {return ac_;}

        // setters
        void setSigma(Real sigma) {sigma_ = sigma;}
        void setAlphaC(Real alphaC) {ac_ = alphaC;}

    protected:
        Real sigma_=0.01;
        Real ac_=0.02;

        // default Dynamic to be false
        bool isDynamic_=false;
        SimulationState& simstate_;
        SimulationBox& simbox_;

        std::string name_;
};

namespace ProbeVolumes
{
    using Base    = ProbeVolume;
    using Key     = std::string;

    using Factory = GenericFactory<Base,Key, ProbeVolumeInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, ProbeVolumeInput&>;
}