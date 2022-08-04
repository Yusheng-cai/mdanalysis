#pragma once

#include "Calculation.h"
#include "LinAlgTools.h"
#include "Bin.h"
#include "SimulationState.h"

#include <vector>
#include <memory>
#include <string>

class CrossSectionalDensity : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;

        // constructor 
        CrossSectionalDensity(const CalculationInput& input);

        // the virtual functions to be overriden 
        virtual void calculate() override;
        virtual void update() override {};
        virtual void finishCalculate() override;

        // function to print histogram
        void printHistogram(std::string name);

    private:
        binptr xbin_;
        binptr ybin_;
        int xbinsize_;
        int ybinsize_;

        std::string resname_;
        std::vector<std::vector<Real>> histogram_;
};