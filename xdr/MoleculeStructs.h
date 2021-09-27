#pragma once
#include "tools/CommonTypes.h"

namespace Molecule
{
    struct atom
    {
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        int residueNumber_;
        std::string residueName_;
        std::string atomName_;

        int atomNumber_;
        Real mass_;
        Real charge_;

        Real3 positions_;
    };

    struct residue
    {
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        std::vector<atom> atoms_;
        std::vector<Real> mass_;
    };
}