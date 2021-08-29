#include "Calculation.h"
#include "CalculationTools.h"
#include "Bin.h"
#include "parallel/OpenMP_buffer.h"
#include "tools/Constants.h"

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
        virtual void printOutput() override;

    private:
        std::string resname1_;
        std::string resname2_;

        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;

        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;

        std::vector<Real> distance_;

        Binptr bins_;

        OpenMP::OpenMP_buffer<std::vector<Real>> distanceBuffer_;

        std::vector<Real> rdf_;

        std::vector<Real> volume_;

        std::string outputName_;
        std::ofstream outputofs_;

        int precision_ = 3;
};