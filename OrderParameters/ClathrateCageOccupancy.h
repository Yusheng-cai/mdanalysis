#include "Calculation.h"
#include "ChillPlus.h"
#include "CellGrid.h"
#include "DensityField.h"
#include "tools/CommonOperations.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <memory>

class ClathrateCageOccupancy : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using densityptr = std::unique_ptr<DensityField>;
        using index3 = CommonTypes::index3;

        ClathrateCageOccupancy(const CalculationInput& input);

        void printAtomsNearby(std::ofstream& ofs);
        void printAtomsWithin(std::ofstream& ofs);

        virtual void calculate() override;
        virtual void finishCalculate() {};
        void update();
        void readFile(const std::string& filename);

    private:
        Real radius_sq_;
        Real radius_;
        Real threshold_;
        std::vector<std::vector<int>> indices_;
        cellptr cell_;

        std::vector<Real> neighbor_water_iter_;
        std::vector<Real> neighbor_water_;

        std::string substrate_name_, index_file_;
        Real3 Ray_{0,0,1};

        densityptr density_;

        int num_within_volume_;
        int num_with_neighbors_;

        Real sigma_ = 0.24;
        Real n_ = 3.0;
        index3 nL_{40,40,40};
        Real isoval_ = 12.5;
        bool pbcMesh_ = true;

        // within methane
        std::vector<int> methane_within_;
        std::vector<int> methane_nearby_;
};