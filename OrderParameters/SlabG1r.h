#pragma once

#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "Bin.h"
#include "LinAlgTools.h"

#include <vector>
#include <string>
#include <memory>

class SlabG1r : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        SlabG1r(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void initializeDistanceCOM();

        void binUsingMinMax();

        void printSlabG1r(std::string name);

    private:
        std::vector<int> distanceCOM_;

        int headindex_ = 1; 
        int tailindex_ = 2;
        std::vector<Real3> uij_;

        std::string residue_;

        Real above_=0;

        int index_=2;

        binptr zbin_, rbin_;
        std::vector<Real> zBinLocation_;
        int numzbins_, numrbins_;

        std::vector<std::vector<int>> IndicesZbin_;

        std::vector<std::vector<Real>> G1r_;
        std::vector<std::vector<Real>> NumberResidueG1r_;

        // are we calculating g1r within bin
        bool g1r_within_zbin_=false;
};