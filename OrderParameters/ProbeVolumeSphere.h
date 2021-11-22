#pragma once
#include "ProbeVolume.h"
#include "SimulationState.h"
#include "IndicatorFunction1d.h"
#include "Calculation.h"

#include <algorithm>
#include <list>
#include <numeric>

class ProbeVolumeSphere: public ProbeVolume
{
    public:
        ProbeVolumeSphere(ProbeVolumeInput& input);
        virtual ~ProbeVolumeSphere(){};

        virtual void setGeometry() override;
        virtual void update() override;

        virtual ProbeVolumeOutput calculate(const Real3& x) const override;

    private:   
        Real sigma_=0.01;
        Real ac_=0.02;

        Real r_;
        Real rmax_;
        Real3 center_;
        IndicatorFunction1d func_;
};