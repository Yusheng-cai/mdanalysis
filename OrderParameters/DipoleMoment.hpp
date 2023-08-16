#include "Calculation.h"
#include "tools/CommonTypes.h"
#include "Bin.h"
#include "SimulationState.h"
#include <string>
#include <cmath>
#include <memory>

class DipoleMoment : public Calculation
{
    public:
        using Binptr = std::unique_ptr<Bin>;

        DipoleMoment(const CalculationInput& input);

        void calculate() override;
        void update() override {}
        void finishCalculate() override;
        
        void printAngleDistribution(std::string name);

    private:
        std::string resname_;
        std::vector<Real3> mu_;

        std::vector<Real> azimuthal_angle_;
        std::vector<Real> zenithal_angle_;

        std::vector<Real> azimuthal_angle_dist_;
        std::vector<Real> zenithal_angle_dist_;

        int numbins_=50;
        Binptr bin_;
        Binptr azi_bin_;
        Range range_ = {{-1, 1}};

        std::vector<int> dipole_index_;
};