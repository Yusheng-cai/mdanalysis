#pragma once
#include "Calculation.h"
#include "tools/Assert.h"
#include "Bin.h"
#include "LinAlgTools.h"
#include "SimulationState.h"

#include <vector>
#include <array>
#include <string>
#include <iomanip>
#include <memory>

class Dipole : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using Range  = std::array<Real,2>;

        Dipole(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
        void printHistogramsquared(std::string name);
        void initializeBins();
        void initializeAtomIndices();

    private:
        std::string residueName_;
        std::vector<Real3> dipoledirection_;

        std::vector<Real> histogram_;
        std::vector<Real> histogramsquared_;

        // The atom indices in which we want to find the dipoles 
        std::vector<int> Atomindices_;
        Real3 direction_ = {{0,0,1}};
        Binptr bin_;
        Binptr binsquared_;

        int numsquaredbin_=20;
};