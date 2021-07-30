#include "OrderParameters.h"
#include "ProbeVolume.h"
#include "parallel/OpenMP_buffer.h"
#include "ProbeVolumeRegistry.h"

class Indus:public OrderParameters
{
    public:
        Indus(const OrderParametersInput& input);
        virtual ~Indus(){};

        virtual void calculate() override;
        virtual void update() override;

        // getters
        Real getN() const {return N_;}
        Real getNtilde() const {return Ntilde_;}
    
    private:
        Real N_;
        Real Ntilde_;

        std::string atomGroupName_;
        std::string pvName_;
        ProbeVolumeRegistry& pv_;
};