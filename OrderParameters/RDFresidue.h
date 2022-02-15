#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "parallel/OpenMP_buffer.h"
#include "tools/Constants.h"
#include "SimulationState.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <array>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>
#include <chrono>

class RDFresidue : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        RDFresidue(const CalculationInput& input);
    
        virtual void calculate() override;
        virtual void finishCalculate() override;

        // print out the radial distribution function that is not normalized by rho
        void printRDFUnnormalized(std::string name);

        // print out the normal radial distribution function
        void printRDF(std::string name);

        // just print the averaged number of neighbors 
        void printNumNeighbors(std::string name);

    private:
        std::string resname1_;
        std::string resname2_;

        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;

        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;

        Binptr bins_;

        std::vector<Real> rdf_;
        Real rho_;

        std::vector<Real> volume_;

        std::string outputName_;
        std::ofstream outputofs_;

        int precision_ = 3;
};