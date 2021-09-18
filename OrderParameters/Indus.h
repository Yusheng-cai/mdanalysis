#pragma once
#include "OrderParameters.h"
#include "ProbeVolume.h"
#include "parallel/OpenMP_buffer.h"

class Indus:public OrderParameters
{
    public:
        using Real3 = CommonTypes::Real3;

        Indus(const OrderParametersInput& input);
        virtual ~Indus(){};

        virtual void calculate() override;
        virtual void update() override;

        // getters
        Real getN() const {return N_;}
        Real getNtilde() const {return Ntilde_;}
        const std::vector<int>& getIndusIndices() const{return indusIndices_;};

        // accessors
        std::vector<int>& accessIndusIndices() {return indusIndices_;} 
    private:
        Real N_;
        Real Ntilde_;

        std::string atomGroupName_;
        std::string pvName_;

        OpenMP::OpenMP_buffer<std::vector<int>> IndusIndicesBuffer_;
        std::vector<int> indusIndices_;

        std::vector<Real3> derivatives_;
};