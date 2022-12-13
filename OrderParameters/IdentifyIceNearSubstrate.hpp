#include "OrderParameters.h"
#include "Calculation.h"
#include "CellGrid.h"

#include <iostream>
#include <memory>
#include <vector>
#include <array>

class IdentifyIceNearSubstrate : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;

        IdentifyIceNearSubstrate(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;
        void readFile(const std::string& name);
        void printNumNeighborWater(std::string name);

    private:
        std::string substrate_name_;
        std::string atom_name_;
        std::string index_file_;

        Real radius_=1.0;
        Real radius_sq_=1.0;

        // set up the cell grid
        cellptr cell_;

        // indices 
        std::vector<std::vector<int>> indices_;

        // number of neighbor waters
        std::vector<Real> neighbor_water_;

        int step_=0;
};