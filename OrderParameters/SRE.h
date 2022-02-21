#pragma once

#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "CellGrid.h"

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
        virtual void finishCalculate() override {};
        virtual void update() override;

        void calculateWithNS();
        void calculateWithoutNS();

        void getNonZeroCharges();

        Real getEnergy() {return energy_;}

        void printEnergyPerIter(std::ofstream& ofs);

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

        Cellptr cell_;

        // indices of the solvent/solute atoms which have nonzero charges
        std::vector<int> NonZeroSolvent_;
        std::vector<int> NonZeroSolute_;

        // mode at which we want to operate at 
        std::string mode_="NS";

        // option to only consider attrative parts 
        bool onlyattrative_=false;
};