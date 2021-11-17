#pragma once
#include "tools/Assert.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "xdr/TopologyReader.h"
#include "xdr/GroFile.h"
#include "AtomGroup.h"
#include "ResidueGroup.h"
#include "tools/CommonTypes.h"
#include "Output_values.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <functional>
#include <memory>

class SimulationState;

struct CalculationInput
{
    ParameterPack& pack_;
    SimulationState& simstate_;
};

class Calculation
{
    public:
        using OutputRegistry = std::map<std::string,OutputValue>; 

        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using outputFunc = std::function<void(std::string)>;
        using perIteroutputFunc = std::function<void(std::ofstream&)>;
        using ofsptr = std::unique_ptr<std::ofstream>;
        using Range = std::array<Real,2>;

        Calculation(const CalculationInput& input);
        virtual ~Calculation(){};

        virtual void update(){};
        virtual void calculate() = 0;
        virtual void finishCalculate() = 0;
        virtual void printOutput();
        virtual void printOutputOnStep();
        virtual void calculateBetaFactors() {};

        void addAtomgroup(std::string name);
        void addResidueGroup(std::string name);

        std::string getName() {return name_;}

        int getNumResidueGroups() const { return ResidueGroups_.size();}
        int getNumAtomGroups() const {return AtomGroups_.size();}

        void registerOutputFunction(std::string name, outputFunc func);
        outputFunc& getOutputByName(std::string name);
        void registerPerIterOutputFunction(std::string name, perIteroutputFunc func);
        perIteroutputFunc& getIterOutputByName(std::string name);

        // register outputs in the output file
        void registerOutputFileOutputs(std::string name, OutputValue::ValueFunction func);
        const OutputRegistry& getOutputRegistry() { return output_;}

        void closeAllOutputPerIter();

        void initializeResidueGroup(const std::string& residueName);

        const ResidueGroup& getResidueGroup(std::string name) const;

        void initializeProbeVolumes();

    protected:
        // output registry 
        OutputRegistry output_;

        std::string name_="c";

        int precision_=3;

        SimulationState& simstate_;

        std::map<std::string,int> MapAtomGroupNameToIndex_;
        std::map<std::string,int> MapResidueGroupNameToIndex_;

        std::vector<AtomGroup*> AtomGroups_; 
        std::vector<ResidueGroup*> ResidueGroups_;

        std::vector<std::string> ResidueNames_;
        std::vector<std::string> AtomGroupNames_;

        ParameterPack& pack_;

        std::vector<int> COMIndices_;
        std::vector<Real3> COM_;

        // map from name to output function
        std::map<std::string, outputFunc> MapNameToOutputFunction_;
        std::map<std::string, perIteroutputFunc> MapNameToPerIterOutput_;

        // vector of outputs
        std::vector<std::string> vectorOutputs_;

        // vector of output file names 
        std::vector<std::string> vectorOutputNames_;

        // vector of ofs for per iter calculation outputs
        std::vector<ofsptr> ofsVector_;

        // vector of per iter outputs
        std::vector<std::string> perIteroutputs_;
        std::vector<std::string> perIteroutputNames_;

        // beta factors 
        std::vector<Real> BetaFactors_;

        // string that holds the names of the probe volumes
        std::vector<std::string> probevolumeNames_;
};

namespace CalculationRegistry
{
    using Key = std::string;
    using Base= Calculation;

    using Factory = GenericFactory<Base, Key, const CalculationInput&>;

    template <typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const CalculationInput&>;
};