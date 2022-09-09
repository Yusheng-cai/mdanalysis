#pragma once

#include "Calculation.h"
#include "CellGrid.h"
#include "SimulationState.h"
#include "tools/Algorithm.h"
#include "tools/CommonOperations.h"
#include "tools/Constants.h"
#include "tools/CommonTypes.h"

#include <memory>
#include <map>
#include <cmath>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <vector>
#include <array>
#include <complex>

enum ChillPlusTypes
{
    Cubic, Hexagonal, Interfacial, Clathrate, Interfacial_Clathrate, Liquid
};

enum BondTypes
{
    eclipse, staggered
};

class ChillPlus : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using INT2    = CommonTypes::index2;

        ChillPlus(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        void printIcetypes(std::string name);

        bool isEclipse(Real cij) {if (cij <= 0.25 && cij >= -0.35) {return true;} return false;}
        bool isStaggered(Real cij) {if (cij <= -0.8) {return true;} return false;}


    private:

        cellptr cell_;
        
        Real solvation_shell_r_=0.35;
        Real solvation_shell_r_squared_;

        int harmonics_degree_=3;
        std::vector<int> m_vec_;
        int num_m_;

        std::string atomgroup_name_;

        std::map<INT2, int> mapBondToIceType_ = 
        {
           {{{0,4}}, ChillPlusTypes::Cubic},
           {{{1,3}}, ChillPlusTypes::Hexagonal},
           {{{0,2}}, ChillPlusTypes::Interfacial},
           {{{1,2}}, ChillPlusTypes::Interfacial}, 
           {{{2,2}}, ChillPlusTypes::Interfacial},
           {{{0,3}}, ChillPlusTypes::Interfacial},
           {{{4,0}}, ChillPlusTypes::Clathrate},
           {{{3,0}}, ChillPlusTypes::Interfacial_Clathrate},
           {{{3,1}}, ChillPlusTypes::Interfacial_Clathrate}
        };

        std::vector<INT2> Bonds_;

        std::vector<std::vector<int>> ice_types_;
};