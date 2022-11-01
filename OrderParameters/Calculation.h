#pragma once
#include "tools/Assert.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "tools/Constants.h"
#include "xdr/TopologyReader.h"
#include "xdr/GroFile.h"
#include "AtomGroup.h"
#include "ResidueGroup.h"
#include "tools/CommonTypes.h"
#include "Output_values.h"
#include "ProbeVolume.h"
#include "tools/CommonOperations.h"
#include "tools/Timer.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <functional>
#include <memory>
#include <math.h>

class SimulationState;

struct CalculationInput{
    ParameterPack& pack_;
    SimulationState& simstate_;
};

class Calculation{
    public:
        using OutputRegistry = std::map<std::string,OutputValue>; 

        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT3 = CommonTypes::index3;
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

        /*
        Function that adds and obtains the atomgroup pointers from the vector
        */
        void addAtomgroup(std::string name);
        const AtomGroup& getAtomGroup(std::string name) const;
        int getNumAtomGroups() const {return AtomGroups_.size();}

        /*
        Function that adds residue group pointer to the local vector 
        */
        void addResidueGroup(std::string name);
        const ResidueGroup& getResidueGroup(std::string name) const;
        int getNumResidueGroups() const { return ResidueGroups_.size();}

        Real3 calcCOM(const Molecule::residue& residues);
        Real3 calcCOM(const Molecule::residue& residues, std::vector<int>& COMIndices);
        Real3 calcCOM(const ResidueGroup& resgroup);

        void registerOutputFunction(std::string name, outputFunc func);
        outputFunc& getOutputByName(std::string name);
        void registerPerIterOutputFunction(std::string name, perIteroutputFunc func);
        perIteroutputFunc& getIterOutputByName(std::string name);

        // register outputs in the output file
        void registerOutputFileOutputs(std::string name, OutputValue::ValueFunction func);
        const OutputRegistry& getOutputRegistry() { return output_;}

        void closeAllOutputPerIter();

        void initializeResidueGroup(const std::string& residueName);
        void initializeResidueGroup(const std::string& residueName, std::string COMName, std::vector<int>& COMIndices, \
        std::vector<Real3>& COM);

        // read any indices , if not read, then the indices will be 0-ResidueLength 
        void ReadResidueIndices(const std::string& residueName, std::string IndicesName, std::vector<int>& Indices);


        // initialize the probe volumes 
        // Probevolumes and not in prove volumes means both have to be satisfied for an atom 
        // to be counted
        void initializeProbeVolumes();
        void initializeNotInProbeVolumes();

        // obtain the indices inside a certain probe volume
        std::vector<int> InsidePVIndices(std::vector<Real3>& pos);
        std::vector<int> InsidePVIndices(std::vector<Real3>& pos, std::vector<int>& outsideIndices);
        bool isInPV(const Real3& pos);
        bool isInPV(const Real3& pos, Real& htildex);

        /*
        Function that gets the name of this particular calculation object 
        */
        std::string getName() {return name_;}

    protected:
        // timer object
        Timer t;
    
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
        std::vector<ProbeVolume*> probevolumes_;

        // string that holds the names of probe volumes to be excluded
        std::vector<std::string> NotInprobevolumeNames_;
        std::vector<ProbeVolume*> NotInprobevolumes_;

        // COM mode 
        std::string COM_mode_="mass";
};

namespace CalculationRegistry
{
    using Key = std::string;
    using Base= Calculation;

    using Factory = GenericFactory<Base, Key, const CalculationInput&>;

    template <typename D>
    using registry_ = RegisterInFactory<Base, D, Key, const CalculationInput&>;
};

namespace CalculationTools
{
    using Real  = CommonTypes::Real;
    using INT3  = std::array<int,3>;
    using Real3 = CommonTypes::Real3;

    Real3 getCOM(const Molecule::residue& residues, const SimulationState& simstate, std::vector<int>& indices);
    Real3 getCOC(const Molecule::residue& residues, const SimulationState& simstate, std::vector<int>& indices);
    Real3 getCOG(const Molecule::residue& residues, const SimulationState& simstate, std::vector<int>& indices);

    // calculate Usr between a pair of residues 
    void CalculateUsrBetweenPair(const Molecule::residue& residue1, const Molecule::residue& residue2, \
                                 const SimulationState& simstate, Real r, \
                                 const std::vector<int>& Usrindices1, const std::vector<int>& Usrindices2, Real sigma, \
                                 Real& attr, Real& repul, Real& total);
    
    // reduce a position to the nearest lattice index
    INT3 NearestLatticeIndex(const Real3& pos, const Real3& dL);
    INT3 NearestLatticeIndex(const Real3& pos, const Real3& dL, const INT3& Lattice_shape);

    // coarse graining function
    Real corase_grain_function(Real r, Real sigma);

    // correct for pbc lattice index
    INT3 correctPBCLatticeIndex(const INT3& latticeIndex, const INT3& lattice_shape);
};