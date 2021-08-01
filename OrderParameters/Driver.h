#include "SimulationState.h"
#include "OrderParameters.h"
#include "ProbeVolume.h"
#include "ProbeVolumeRegistry.h"
#include "xdr/XdrWrapper.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "AtomGroup.h"
#include "tools/Registry.h"
#include "Output_values.h"
#include "Output_files.h"

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <chrono>

class Driver
{
    public:
        using ProbeVolumePtr = std::unique_ptr<ProbeVolume>;
        using OPptr          = std::unique_ptr<OrderParameters>;
        using XdrPtr         = std::unique_ptr<XdrWrapper>;
        using Real           = CommonTypes::Real;
        using Real3          = CommonTypes::Real3;
        using VectorReal3    = CommonTypes::VectorReal3;
        using outputptr      = std::unique_ptr<OutputStream>;

        Driver(std::string filename);
        ~Driver(){};

        void initializeXdr(const ParameterPack*);
        void initializeOP(const std::vector<const ParameterPack*>&);
        void initializeProbeVolume(const std::vector<const ParameterPack*>&);
        void initializeAtomGroups(const std::vector<const ParameterPack*>&);
        void intializeOutputFiles(const std::vector<const ParameterPack*>&);

        void RegisterOuputValues();
        const OutputValue& getOutputValue(std::string name) const;

        Real getTime() const {return Xdr_->getTime();}
        int getStep() const {return Xdr_->getStep();}
        int getNframes() const {return Xdr_->getNframes();}

        void update();
        void calculate();

        bool isActive(){return is_Active_;}

    private:
        ParameterPack pack_;
        std::map<std::string, ProbeVolumePtr> MapName2PV_;
        std::vector<OPptr> OP_;
        XdrPtr Xdr_;
        SimulationState simstate_;
        ProbeVolumeRegistry pv_registry_;

        VectorReal3 total_atom_positions_;

        std::vector<std::string> VectorAgNames_;

        Registry<std::string, OutputValue> outputValueRegistry_;

        std::vector<outputptr> OutputFiles_;

        bool is_Active_ = true;
};