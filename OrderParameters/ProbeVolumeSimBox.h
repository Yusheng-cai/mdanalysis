#pragma once
#include "ProbeVolume.h"
#include "tools/CommonTypes.h"

#include <array>
#include <vector>

class ProbeVolumeSimBox:public ProbeVolume
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Real2= CommonTypes::Real2;

        ProbeVolumeSimBox(ProbeVolumeInput& input);
        virtual ~ProbeVolumeSimBox(){};

        // This calculates per atom basis of indicator function
        virtual ProbeVolumeOutput calculate(const Real3& x) const override; 
        virtual void setGeometry() override {};


    private:
};