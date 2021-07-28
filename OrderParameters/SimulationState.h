#pragma once
#include "tools/CommonTypes.h"
#include "SimulationBox.h"

// This is a class that keeps track of the simulation progression
class SimulationState
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Matrix=CommonTypes::Matrix;

        SimulationState(){};

        void setBox(Matrix boxMat){box_.setBoxMatirx(boxMat);}
        void setTime(Real time){time_ = time;}
        void setStep(int step){step_ = step;}

        // getters
        SimulationBox& getSimulationBox(){return box_;}
        const SimulationBox& getSimulationBox() const{return box_;}
        Real getTime() const{return time_;}
        int getStep() const {return step_;}

    private:
        Real time_;
        int step_;

        SimulationBox box_;
};