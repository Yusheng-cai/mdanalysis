#pragma once
#include "tools/CommonTypes.h"
#include "tools/Constants.h"

#include <cmath>

class IndicatorFunction
{
    public:
        using Real = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;

        IndicatorFunction(){};
        IndicatorFunction(Real sigma_, Real ac_);
        virtual ~IndicatorFunction(){};

        // h_a of a single dimension in INDUS indicator function
        virtual void calculate(const Real& x, Real& h_x, Real& htilde_x, Real& dhtilde_dx) const = 0;
        virtual void setLimits() = 0;
        virtual void calculateFactors();

        // truncated Gaussian phi function, this is a zero centered truncated Gaussian 
        Real truncatedGaussian(Real alpha) const;

        // getters
        Real getSigma() const{return sigma_;}
        Real getAlphaC() const{return ac_;}
        
        // setters
        void setSigma(Real sigma) {sigma_ = sigma;}
        void setAlphaC(Real ac) {ac_ = ac;}
    
    protected:
        Real sigma_;
        Real sigma2_;
        Real ac_;
        Real ac2_;

        Real k_;
        Real k1_;
        Real k2_;
};