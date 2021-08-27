#include "Calculation.h"
#include "CalculationTools.h"
#include "Bin.h"
#include "parallel/OpenMP_buffer.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <array>
#include <string>
#include <cmath>
#include <memory>

class RDFresidue : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;
        RDFresidue(const CalculationInput& input);
    
        virtual void calculate() override;
        virtual void finishCalculate() override {};

    private:
        std::string resname1_;
        std::string resname2_;

        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;

        bool COMIndices1Read_=false;
        bool COMIndices2Read_=false;

        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;

        std::vector<Real> distance_;

        Binptr bins_;

        OpenMP::OpenMP_buffer<std::vector<Real>> distanceBuffer_;

        std::vector<int> numCountsPerBin_;
};