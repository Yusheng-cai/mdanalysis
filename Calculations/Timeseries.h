#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <numeric>

class Timeseries
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        Timeseries(const ParameterPack& pack);
    private:
        std::vector<Real> data_;
}