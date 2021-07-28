#include "tools/CommonTypes.h"

#include <cmath>

class IndicatorFunction
{
    public:
        using Real = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;

        IndicatorFunction();
        virtual ~IndicatorFunction(){};
        virtual void calculate(const Real3& x, Real& h_x, Real& htilde_x, Real3& dhtilde_dx) = 0;

        // getters
        Real getSigma() const{return sigma_;}
        Real getAlphaC() const{return ac_;}
        
        // setters
        void setSigma(Real sigma) {sigma_ = sigma;}
        void setAlphaC(Real ac) {ac_ = ac;}
    
    protected:
        Real sigma_;
        Real ac_;
};