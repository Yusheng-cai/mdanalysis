#pragma once

#include "ProbeVolume.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "ProbeVolumeCylinder.h"
#include "AtomGroup.h"
#include "LinAlgTools.h"

#include <memory>
#include <vector>
#include <array>
#include <string>

class ProbeVolumeTiltedCylinder : public ProbeVolume
{
    public:
        using cylinderPtr = std::unique_ptr<ProbeVolumeCylinder>; 
        using Matrix      = CommonTypes::Matrix;
        ProbeVolumeTiltedCylinder(ProbeVolumeInput& input);

        virtual ProbeVolumeOutput calculate(const Real3& x) const;
        virtual void update() override;
        virtual void setGeometry() override;
        

    private:
        cylinderPtr cylinder_;

        Real3 Zvector_;

        // define base and top
        Real3 base_;
        Real3 top_;

        // distance 
        Real distance_;

        // reference vector 
        Real3 refVector = {{1,0,0}};

        // The rotation matrix 
        Matrix rotationMat_;

        Real radius_;
};