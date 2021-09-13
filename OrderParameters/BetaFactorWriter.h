#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"

#include <iostream>
#include <string>
#include <vector>

class BetaFactorWriter
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        BetaFactorWriter(ParameterPack& param);

        void write(int frameNum, const std::vector<Real>& data);
    
    private:
        std::ofstream ofs_;
        std::string filename_;
};