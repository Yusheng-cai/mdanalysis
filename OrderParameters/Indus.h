#include "OrderParameters.h"
#include "ProbeVolume.h"

class Indus:public OrderParameters
{
    public:
        Indus(const OrderParametersInput& input);
        virtual ~Indus(){};

        virtual void calculate() override{};
    
    private:
        Real N_;
        Real Ntilde_;
};