#include "ProbeVolume.h"
#include "IndicatorFunction1d.h"

class ProbeVolumeSphere: public ProbeVolume
{
    public:
        ProbeVolumeSphere(ProbeVolumeInput& input);
        virtual ~ProbeVolumeSphere(){};

        virtual void setGeometry() override;

        virtual ProbeVolumeOutput calculate(const Real3& x) const override;

    private:   
        Real sigma_=0.01;
        Real ac_=0.02;

        Real r_;
        Real rmax_;
        Real3 center_;
        IndicatorFunction1d func_;
};