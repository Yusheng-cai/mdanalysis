#pragma once

#include "Calculation.h"
#include "LinAlgTools.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"

#include <vector>
#include <array>
#include <map>
#include <string>

// Deuterium order parameter is defined as 
// Sch = 2/3 * Sxx + 1/3 * Syy
// Sxx = 1/2 \langle 3*cos^{2}(\theta_{x}) - 1 \rangle

class DeuteriumOP : public Calculation
{
    public:
        using Matrix = CommonTypes::Matrix;
        DeuteriumOP(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() {};
        virtual void finishCalculate() override;

        void printDeuteriumOP(std::string name);
    
    private:    
        std::string residue_;
        std::vector<int> CarbonIndices_;
        int numCarbons_;

        std::vector<Real> avgSCD_;
        std::vector<Real> avgSCC_;
        std::vector<Real> avgSZZ_;

        Real3 SurfaceNormal_={{0,0,1}};
};