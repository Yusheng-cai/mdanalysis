#pragma once 

#include "Calculation.h"
#include "Bin.h"
#include "tools/Constants.h"
#include "tools/CommonTypes.h"

#include <memory>

class RDF2d : public Calculation
{
    public:
        using INT2   = CommonTypes::index2;
        using Binptr = std::unique_ptr<Bin>;

        RDF2d(const CalculationInput& input);

        virtual void calculate() override;

    private:
        std::string resname1_;
        std::string resname2_;
        Binptr bins_;
        std::vector<Real> Area_;

        std::vector<int> COMIndices1_;
        std::vector<int> COMIndices2_;

        std::vector<Real> rdf_;

        std::vector<Real3> COM1_;
        std::vector<Real3> COM2_;

        INT2 position_index_;
};