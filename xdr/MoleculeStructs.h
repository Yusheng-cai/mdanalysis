#pragma once
#include "tools/CommonTypes.h"

namespace Molecule
{
    struct atom
    {
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        int resnum_;
        std::string resname_;
        std::string atomName_;

        int atomNumber_;
        Real mass_;
        Real charge_;

        std::string type_;
        Real3 positions_;
    };

    struct residue
    {
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        std::vector<atom> atoms_;
    };

    struct AtomType
    {
        using Real = CommonTypes::Real;

        std::string type_;
        std::string resname_;
        std::string atomName_;
        int atomNumber_;
        Real charge_;
        Real mass_;
    };
}