#pragma once
#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"
#include "Bin.h"
#include "algorithm"

#include <memory>
#include <string>
#include <map>


class SlabDimer : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;

        SlabDimer(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override {};
        virtual void finishCalculate() override;
        void binUsingMinMax();

        void printDimerRatio(std::string name);
        void printDimerOrientation(std::string name);
        void printMonomerOrientation(std::string name);

    private:
        Binptr zBin_, tBin_;
        bool usingMinMax_=false;
        int  numzbins_=30, numtbins_=30;
        Real above_=-1e20;
        Real below_=1e20;

        std::vector<std::vector<int>> BinIndices_;
        std::vector<std::vector<int>> DimerIndices_, MonomerIndices_;
        std::vector<int> numCOM_;

        std::string resname_;
        int headindex_=1, tailindex_=2;

        Real alignment_cutoff_, distance_cutoff_, distance_cutoff_B1B2_;
        std::string direction_;
        int index_;

        // Map direction string to index
        std::map<std::string, int> MapNameToDirection = 
        {
            {"x", 0},
            {"y", 1},
            {"z", 2}
        };

        std::vector<Real> BinLocation_;

        std::vector<Real3> COMB1_ , COMB2_;
        std::vector<Real3> uij_;
        std::vector<int> num_dimers_;

        std::vector<int> COMIndicesB1_={3,4,5,6,7,8};
        std::vector<int> COMIndicesB2_={9,10,11,12,13,14};

        std::vector<std::vector<Real>> orientation_dimer_;
        std::vector<std::vector<Real>> orientation_monomer_;

        std::vector<Real> ratio_dimer_;
};
