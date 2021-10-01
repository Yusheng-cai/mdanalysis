#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"
#include "Pcost.h"

#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>

// calculate distribution of cosine theta of the angles in a probe volume

class Pcost2 : public Pcost
{
    public:
        Pcost2(const CalculationInput& input);

        virtual void calculate() override;
};