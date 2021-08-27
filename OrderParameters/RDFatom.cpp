#include "RDFatom.h"

RDFatom::RDFatom(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("atomgroup1", ParameterPack::KeyType::Optional,atomgroup1Name_);
    input.pack_.ReadString("atomgroup2", ParameterPack::KeyType::Optional,atomgroup2Name_);
}

void RDFatom::calculate()
{

}