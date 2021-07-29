#include "OrderParameters.h"
#include "ProbeVolume.h"

class Indus:public OrderParameters
{
    public:
        Indus(const OrderParametersInput& input);
        virtual ~Indus(){};

        virtual void calculate() override{};

        Real getN() const {return N_;}
        Real getNtilde() const {return Ntilde_;}
    
    private:
        Real N_;
        Real Ntilde_;
};