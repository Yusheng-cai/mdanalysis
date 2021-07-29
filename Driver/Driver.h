#include "OrderParameters/SimulationState.h"
#include "OrderParameters/OrderParameters.h"
#include "OrderParameters/ProbeVolume.h"
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
        using statePtr = std::unique_ptr<SimulationState>;

        Driver(std::string filename);
        ~Driver(){};

        void initializeXdr(const ParameterPack*);
        void initializeOP(std::vector<const ParameterPack*>);
        void initializeProbeVolume(std::vector<const ParameterPack*>);
    private:
        ParameterPack pack_;
        std::vector<ProbeVolumePtr> PV_;
        std::vector<OPptr> OP_;
        XdrPtr Xdr_;
        statePtr SimState_;
};