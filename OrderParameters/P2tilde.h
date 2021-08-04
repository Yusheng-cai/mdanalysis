#include "OrderParameters.h"
#include "Qtensor.h"
#include "tools/CommonTypes.h"
#include "ProbeVolumeRegistry.h"

#include <vector>
#include <string>

class P2tilde: public OrderParameters
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real   = CommonTypes::Real;
        using Real3  = CommonTypes::Real3;

        P2tilde(const OrderParametersInput& input);
        virtual ~P2tilde(){};

        virtual void calculate() override;
        virtual void update() override;

        Real getP2tilde() const {return p2tilde_;}
        Real getQxx() const {return Qtensor_[0][0];}
        Real getQxy() const {return Qtensor_[0][1];}
        Real getQxz() const {return Qtensor_[0][2];}
        Real getQyy() const {return Qtensor_[1][1];}
        Real getQyz() const {return Qtensor_[1][2];}

    private:
        Matrix Qtensor_;

        Real p2tilde_;

        // obtain the registry for probevolume
        ProbeVolumeRegistry& pv_;

        // obtain the name of the probevolume
        std::string pvname_;

        // obtain the name of the head and tail atom groups
        std::string headgroupname_;
        std::string tailgroupname_;

        // the indices of the atoms within the probeVolume
        std::vector<int> IndusIndices_;
        OpenMP::OpenMP_buffer<std::vector<int>> IndusIndicesbuffer_;

        // Ntilde and N of the system
        Real Ntilde_, N_;
};