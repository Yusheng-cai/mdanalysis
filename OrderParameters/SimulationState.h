#pragma once
#include "tools/CommonTypes.h"
#include "SimulationBox.h"
#include "AtomGroup.h"
#include "ResidueGroup.h"
#include "ProbeVolume.h"
#include "Calculation.h"

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
        using ProbeVolumePtr = std::unique_ptr<ProbeVolume>;

        SimulationState(){};
        ~SimulationState(){};

        void setSimulationBox(Matrix boxMat){box_.setBoxMatrix(boxMat);}
        void setTime(Real time){time_ = time;}
        void setStep(int step){step_ = step;}
        void setTotalFrames(int totalFrames){totalframes_ = totalFrames;}
        void setTotalNumberAtoms(int totalAtoms){totalNumberAtoms_=totalAtoms;}
        void setFrameNumber(int FrameNum) {FrameNumber_ = FrameNum;}
        void setTotalAtomPos(const std::vector<Real3>& totalAtomPos) 
        {totalAtomPos_.clear();
         totalAtomPos_.insert(totalAtomPos_.end(), totalAtomPos.begin(), totalAtomPos.end());}

        // getters
        SimulationBox& getSimulationBox(){return box_;}
        const SimulationBox& getSimulationBox() const{return box_;}
        const std::vector<Real3>& getTotalAtomPos() const {return totalAtomPos_;}
        int getTotalFrames() const {return totalframes_;}
        Real getTime() const{return time_;}
        int getStep() const {return step_;}
        int getTotalNumberAtoms() const {return totalNumberAtoms_;}
        int getFrameNumber() const {return FrameNumber_;}

        // registers AtomGroups
        void registerAtomGroup(std::string name, AtomGroup& ag);

        // registers ResidueGroups
        void registerResidueGroup(std::string name, ResidueGroup& res);

        // registers Calculation
        void registerCalculation(std::string name, const Calculation& calc);

        // get calculations 
        const Calculation& getCalculation(std::string name);

        // get AtomGroup reference by name
        AtomGroup& getAtomGroup(std::string name);
        const std::map<std::string, AtomGroup>& getAtomGroupRegistry() const{ return MapName2AtomGroup_;}

        // get ResidueGroups by name
        ResidueGroup& getResidueGroup(std::string name);
        const std::map<std::string, ResidueGroup>& getResidueGroupRegistry() const {return MapName2ResidueGroup_;}

        // register probe volume
        void registerProbeVolume(const std::string name, const ProbeVolume* pv_ptr);
        ProbeVolume& getProbeVolume(const std::string name);

        ProbeVolume& accessProbeVolume(const std::string name);

        // obtain the name to residue map
        const std::map<std::string, ResidueGroup>& getResidueMap() const {return MapName2ResidueGroup_;}

        // obtain the name to atom group map
        const std::map<std::string, AtomGroup>& getAtomGroupMap() const {return MapName2AtomGroup_;}
        

    private:
        Real time_;
        int step_;

        SimulationBox box_;
        std::map<std::string,AtomGroup> MapName2AtomGroup_;
        std::map<std::string,ResidueGroup> MapName2ResidueGroup_;
        std::map<std::string, ProbeVolumePtr> MapName2PV;  
        std::map<std::string, Calculation*> MapName2Calculation_;

        int totalframes_;

        int totalNumberAtoms_;

        int FrameNumber_;

        std::vector<Real3> totalAtomPos_;
};