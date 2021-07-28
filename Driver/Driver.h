#include "OrderParameters/SimulationState.h"
#include "OrderParameters/OrderParameters.h"
#include "OrderParameters/ProbeVolume.h"
#include "xdr/XdrWrapper.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "OrderParameters/SimulationState.h"

#include <string>
#include <memory>
#include <vector>

class Driver
{
    public:
        using ProbeVolumePtr = std::unique_ptr<ProbeVolume>;
        using OPptr = std::unique_ptr<OrderParameters>;
        using XdrPtr = std::unique_ptr<XdrWrapper>;
        using statePtr = std::auto_ptr<SimulationState>;

        Driver(std::string filename);
        ~Driver(){};

        void initializeXdr(std::vector<const ParameterPack*>);
        void initializeOP(std::vector<const ParameterPack*>);
        void initializeProbeVolume(std::vector<const ParameterPack*>);
    private:
        ParameterPack pack;
        std::vector<ProbeVolumePtr> PV;
        std::vector<OPptr> OP;
        std::vector<XdrPtr> Xdr;
        std::vector<statePtr> sim_state_;
};