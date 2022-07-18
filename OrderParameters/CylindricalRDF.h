#pragma once 

#include "Calculation.h"
#include "Bin.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include "tools/Constants.h"
#include "tools/CommonTypes.h"

#include <memory>

class CylindricalRDF : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        using Matrix = CommonTypes::Matrix;

        CylindricalRDF(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printCylindricalRDF(std::string name);
        void printCylindricalUsr(std::string name);
        void printCylindricalUsrAttr(std::string name);

    private:
        // the h direction or z direction is always [0,0,1]
        Real3 direction_={{0,0,1}};

        // residue names
        std::string residue1_;
        std::string residue2_;

        // COMindices 
        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;

        // COM 
        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;

        // rbin and zbin
        binptr rbin_;
        binptr zbin_;
        int numrbin_;
        int numzbin_;

        // number of residues 
        int numres_;

        // average rho
        Real rho_;
        bool interface_=false;

        // use director as array
        bool usedirector_=false;
        int headindex_ = 1;
        int tailindex_ = 2;

        // cylindrical rdf
        std::vector<std::vector<Real>> cylindrical_rdf_;

        // cylindrical volume 
        std::vector<std::vector<Real>> cylindrical_volume_;

        // cylindrical property <\delta(r-rij)\delta(h-hij) property(rij)>
        std::vector<std::vector<Real>> cylindrical_Usr_;
        std::vector<std::vector<Real>> cylindrical_Usr_attr_;
        std::vector<std::vector<Real>> cylindrical_Usr_repul_;
        Real Usr_cutoff_=1.2;
        bool calcUsr_=false;
        Real beta_;

        // Usr stuff
        std::vector<int> UsrIndices1_;
        std::vector<int> UsrIndices2_;
};