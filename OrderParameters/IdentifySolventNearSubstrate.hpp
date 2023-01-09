#include "Calculation.h"
#include "ProbeVolumeSphere.h"
#include "tools/CommonTypes.h"
#include "CellGrid.h"
#include "tools/CommonOperations.h"

#include <string>
#include <memory>

class IdentifySolventNearSubstrate : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using sphereptr = std::unique_ptr<ProbeVolumeSphere>;

        IdentifySolventNearSubstrate(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        void printPerIterAtomsInSphere(std::ofstream& ofs);

    private:
        std::string solvent_name_;
        std::string substrate_name_;

        int num_substrate_atoms_;
        std::vector<Real> num_atoms_in_sphere_;
        std::vector<Real> PerIter_num_atoms_in_sphere_;

        Real radius_;
        Real radius_sq_;

        cellptr cell_;
};