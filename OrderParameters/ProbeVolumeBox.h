#pragma once
#include "ProbeVolume.h"
#include "tools/CommonTypes.h"
#include "IndicatorFunction2d.h"
#include "SimulationState.h"

#include <array>
#include <vector>

class ProbeVolumeBox:public ProbeVolume
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Real2= CommonTypes::Real2;

        ProbeVolumeBox(ProbeVolumeInput& input);
        virtual ~ProbeVolumeBox(){};

        // This calculates per atom basis of indicator function
        virtual ProbeVolumeOutput calculate(const Real3& x) override; 
        void setGeometry();

    private:
        Real sigma_=0.01;
        Real ac_=0.02;

        Real2 xrange_;
        Real2 yrange_;
        Real2 zrange_;
     
        Real dx_,dy_,dz_;
        Real3 center_;
        std::array<IndicatorFunction2d,3> func_;

        SimulationState& state;
};