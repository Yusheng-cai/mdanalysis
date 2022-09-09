#pragma once

#include "Calculation.h"
#include "SimulationState.h"
#include "tools/CommonOperations.h"
#include "LinAlgTools.h"

#include <vector>
#include <array>
#include <string>

class SmecticOrder : public Calculation
{
    public:
        SmecticOrder(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;

    private:
        std::string residuegroup_;

        int headindex_=1, tailindex_=2;
        std::vector<Real3> uij_;
};