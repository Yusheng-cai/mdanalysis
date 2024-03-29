#pragma once

#include "SimulationState.h"
#include "OrderParameters.h"
#include "ProbeVolume.h"
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
#include "Calculation.h"
#include "ResidueGroup.h"

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
        using calcptr        = std::unique_ptr<Calculation>;
        using XdrPtr         = std::unique_ptr<XdrWrapper>;
        using Real           = CommonTypes::Real;
        using Real3          = CommonTypes::Real3;
        using VectorReal3    = CommonTypes::VectorReal3;
        using outputptr      = std::unique_ptr<OutputStream>;

        Driver(std::string filename, CommandLineArguments& cmd);
        ~Driver(){};

        // initialize the xdr file
        void initializeXdr();

        // initialize the topology file
        void initializeTop();

        // initialize OrderParameters
        void initializeOP();

        // initialize the calculations
        void initializeCalculation();
        void initializeProbeVolume();
        void initializeAtomGroups();
        void initializeOutputFiles();
        void initializeGroFile();
        void initializeDriverPack();
        void initializeResidueGroups();

        void RegisterOuputValues();
        const OutputValue& getOutputValue(std::string name) const;

        // function that checks whether or not a step is used to be calculated
        bool isValidStep(int step);

        Real getTime() const {return Xdr_->getTime();}
        int getStep() const {return Xdr_->getStep();}
        int getNframes() const {return Xdr_->getNframes();}

        bool readFrame(int FrameNum);
        void update();
        void calculate();
        void finishCalculate();
        void printOutput();

        bool isActive(){return is_Active_;}

    private:
        ParameterPack pack_;

        // probe volume names 
        std::vector<std::string> ProbeVolumeNames_;

        // Vector of OrderParameters
        std::vector<OPptr> OP_;

        // Vector of calculation objects
        std::vector<calcptr> Calc_;

        // Xdr file pointer
        XdrPtr Xdr_;

        // simulation state object
        SimulationState simstate_;

        VectorReal3 total_atom_positions_;

        std::vector<std::string> VectorAgNames_;
        std::vector<std::string> VectorResNames_;

        Registry<std::string, OutputValue> outputValueRegistry_;

        std::vector<outputptr> OutputFiles_;

        GroFile grofile_;

        bool is_Active_ = true;

        // absolute path of the current working directory
        std::string apath_;

        int startingFrame_ = 1;
        int skip_ = 0;

        // Topolgy obj
        TopologyReader top_;
        bool topology_read_=false;
};