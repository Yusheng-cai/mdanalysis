#pragma once 

#include "Calculation.h"
#include "Bin.h"
#include "SimulationState.h"

#include <vector>
#include <array>
#include <map>
#include <memory>

class EndToEndDistance : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;

        EndToEndDistance(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        void printHistogram(std::string name);
    
    private:
        int head_index_;
        int tail_index_;

        std::vector<Real> histogram_;

        std::string resname_;

        binptr bin_;
        int numbins_;
};