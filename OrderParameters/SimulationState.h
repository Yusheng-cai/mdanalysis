#pragma once
#include "tools/CommonTypes.h"
#include "SimulationBox.h"
#include "AtomGroup.h"

#include <string>
#include <map>
#include <memory>

// This is a class that keeps track of the simulation progression
class SimulationState
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Matrix=CommonTypes::Matrix;

        SimulationState(){};
        ~SimulationState(){};

        void setBox(Matrix boxMat){box_.setBoxMatrix(boxMat);}
        void setTime(Real time){time_ = time;}
        void setStep(int step){step_ = step;}

        // getters
        SimulationBox& getSimulationBox(){return box_;}
        const SimulationBox& getSimulationBox() const{return box_;}
        Real getTime() const{return time_;}
        int getStep() const {return step_;}

        // registers AtomGroups
        void registerAtomGroup(std::string name, AtomGroup& ag);

        // get AtomGroup reference by name
        const AtomGroup& getAtomGroup(std::string name) const;

    private:
        Real time_;
        int step_;

        SimulationBox box_;
        std::map<std::string, AtomGroup> MapName2AtomGroup_;
};