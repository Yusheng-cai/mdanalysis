#include "SimulationState.h"
#include "OrderParameters.h"
#include "ProbeVolume.h"
#include "ProbeVolumeRegistry.h"
#include "xdr/XdrWrapper.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <string>
#include <memory>
#include <vector>
#include <map>

class Driver
{
    public:
        using ProbeVolumePtr = std::unique_ptr<ProbeVolume>;
        using OPptr = std::unique_ptr<OrderParameters>;
        using XdrPtr = std::unique_ptr<XdrWrapper>;
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        Driver(std::string filename);
        ~Driver(){};

        void initializeXdr(const ParameterPack*);
        void initializeOP(std::vector<const ParameterPack*>&);
        void initializeProbeVolume(std::vector<const ParameterPack*>&);

        Real getTime() const {return Xdr_->getTime();}
        int getStep() const {return Xdr_->getStep();}

    private:
        ParameterPack pack_;
        std::map<std::string, ProbeVolumePtr> MapName2PV_;
        std::vector<OPptr> OP_;
        XdrPtr Xdr_;
        SimulationState simstate_;
        ProbeVolumeRegistry pv_registry_;
};