#include "Indus.h"

namespace OrderParametersRegistry
{
    registry_<Indus> register_Indus("indus");
}

Indus::Indus(const OrderParametersInput& input)
:OrderParameters(input)
{
    std::string pv_name;
    input.pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pv_name);

    auto& pv = input.pv_registry_.getProbeVolume(pv_name);
    pv.setGeometry();

    registerOutput("n", [this](void)->Real {return this->getN();});
    registerOutput("ntilde", [this](void)->Real {return this->getNtilde();});
}