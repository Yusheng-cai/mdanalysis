#pragma once

#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "CellGrid.h"
#include "tools/InputParser.h"

#include <vector>
#include <memory>
#include <cmath>
#include <string>
#include <array>

class SRE : public Calculation
{
    public:
        using index3 = CommonTypes::index3;
        using Cellptr= std::unique_ptr<CellGrid>;

        SRE(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;
        virtual void update() override;

        void calculateWithNS();
        void calculateWithoutNS();

        void initializeSoluteSolventIndices();

        void getInsidePVIndices();

        const std::vector<Real> getPerAtomContributionIter() const {return PerAtomContributionIter_;}

        Real getEnergy() {return energy_;}
        Real getRepulsiveEnergy() {return repulsive_energy_;}
        Real getAttractiveEnergy() {return attractive_energy_;}

        void printEnergyPerIter(std::ofstream& ofs);
        void printEnergyPerAtom(std::string name);

    private:
        std::string SolventName_;
        std::string SoluteName_;
        Real epsilon_;
        Real alpha_;
        Real cutoff_=1.2;
        Real cutoffsq_;

        // 1/(4 pi epsilon) as defined by gromac --> kJ nm /(mol e2) 
        Real factor_ = 138.935458;

        Real energy_=0.0;
        Real repulsive_energy_=0.0;
        Real attractive_energy_=0.0;

        Cellptr cell_;

        bool self_interaction_=false;

        // indices of the solvent/solute atoms which have nonzero charges
        std::vector<int> NonZeroSolvent_;
        std::vector<int> NonZeroSolute_;

        // mode at which we want to operate at 
        std::string mode_="NS";

        // per atom contribution to the SRE energy 
        std::vector<Real> PerAtomContribution_;
        std::vector<Real> PerAtomContributionIter_;
        std::vector<Real> PerResidueContribution_;

        // solute indices per residue 
        std::vector<int> SoluteIndices_;
        std::vector<int> SolventIndices_;

        // indus indicator
        std::vector<Real> nonzero_indus_indicators_;
};