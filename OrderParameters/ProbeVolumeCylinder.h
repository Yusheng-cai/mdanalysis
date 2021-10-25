#pragma once
#include "ProbeVolume.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "IndicatorFunction1d.h"
#include "IndicatorFunction2d.h"

#include <vector>
#include <array>

// This is a cylinder class that is only in the z direction, this is for ease of use for tilted cylinder
// This cylinder is centered at [0,0,0]
class ProbeVolumeCylinder : public ProbeVolume
{
    public:
        using Range = CommonTypes::Real2;

        
        ProbeVolumeCylinder(ProbeVolumeInput& input);

        // The x here is already pbc corrected and shifted to the correct frame of reference 
        virtual ProbeVolumeOutput calculate(const Real3& x) const;

        void setGeometry(Real Rmax, Real zmax, Real ac, Real sigma);
    
    private:
        IndicatorFunction2d zfunc_;
        IndicatorFunction1d rfunc_;
        
        Range rrange_;
        Range zrange_;

        // The center of the cylinder, all the atoms will be shifted accordingly
        Real3 center_ = {{0,0,0}};
};