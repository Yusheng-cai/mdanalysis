#pragma once

#include "Calculation.h"
#include "SimulationState.h"
#include "Bin.h"

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <map>
#include <functional>

class SlabAtomicProperty : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        using calcfunc = std::function<std::vector<Real>(std::vector<Real3>&)>;

        SlabAtomicProperty(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;
        std::vector<Real> CalculateDensity(std::vector<Real3>& COM);

        void registerCalcFunc(std::string name, calcfunc func);

        void printProperty(std::string name);
    
    private:
        std::string residueName_;
        std::string property_;
        int direction_=2;

        // bin information
        Binptr bin_;
        int numbins_;
        Real dz_;

        // map from name to calculation function
        std::map<std::string, calcfunc> MapNameToCalculationFunc_;

        // quantity as a function of z
        std::vector<Real> AtomicProperty_;

};