#pragma once

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
#include "xdr/GroFile.h"
#include "tools/CommandLineArguments.h"
#include "tools/FileSystem.h"

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

        Driver(std::string filename, CommandLineArguments& cmd);
        ~Driver(){};

        void initializeXdr(const ParameterPack*);
        void initializeOP(const std::vector<const ParameterPack*>&);
        void initializeProbeVolume(const std::vector<const ParameterPack*>&);
        void initializeAtomGroups(const std::vector<const ParameterPack*>&);
        void initializeOutputFiles(const std::vector<const ParameterPack*>&);
        void initializeGroFile(const ParameterPack*);

        void RegisterOuputValues();
        const OutputValue& getOutputValue(std::string name) const;

        Real getTime() const {return Xdr_->getTime();}
        int getStep() const {return Xdr_->getStep();}
        int getNframes() const {return Xdr_->getNframes();}

        bool readNextFrame();
        void update();
        void calculate();

        bool isActive(){return is_Active_;}

    private:
        ParameterPack pack_;
        std::map<std::string, ProbeVolumePtr> MapName2PV_;

        // Vector of Order Parameters
        std::vector<OPptr> OP_;

        // Xdr file pointer
        XdrPtr Xdr_;

        // simulation state object
        SimulationState simstate_;

        // probe Volume registry that keeps track of all the probe volumes
        ProbeVolumeRegistry pv_registry_;

        VectorReal3 total_atom_positions_;

        std::vector<std::string> VectorAgNames_;

        Registry<std::string, OutputValue> outputValueRegistry_;

        std::vector<outputptr> OutputFiles_;

        GroFile grofile_;

        bool is_Active_ = true;

        // absolute path of the current working directory
        std::string apath_;
};