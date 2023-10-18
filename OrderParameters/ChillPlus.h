#pragma once

#include "Calculation.h"
#include "CellGrid.h"
#include "SimulationState.h"
#include "tools/Algorithm.h"
#include "tools/CommonOperations.h"
#include "tools/Constants.h"
#include "tools/CommonTypes.h"
#include "DensityField.h"

#include <memory>
#include <chrono>
#include <map>
#include <cmath>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <vector>
#include <array>
#include <complex>

enum ChillPlusTypes{
    Cubic, Hexagonal, Interfacial, Clathrate, Interfacial_Clathrate, Surface, Liquid
};

enum BondTypes{
    eclipse, staggered
};

class ChillPlus : public Calculation{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using INT2    = CommonTypes::index2;
        using densityptr = std::unique_ptr<DensityField>;

        ChillPlus(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        Real getNumIceLikeAtoms() const {return num_ice_like_atoms_;}
        Real getNumIceLikeChillPlus() const {return num_ice_like_atoms_before_correction_;}

        void CorrectIceLikeAtomsBasedOnSurface();
        void CorrectForTrueIce();
        void CorrectForTrueClathrate();

        void ShiftTriangleWithRef(Real3& A, Real3& B, Real3& C, Real3& ref);

        void printIceLikeAtoms(std::ofstream& ofs);
        void printClathrateLikeAtoms(std::ofstream& ofs);
        void printIceTypeNum(std::ofstream& ofs);
        void printAddedIce(std::ofstream& ofs);

        bool isEclipse(Real cij) {if (cij <= 0.18 && cij >= -0.45) {return true;} return false;}
        bool isStaggered(Real cij) {if (cij <= -0.8 && cij >= -1) {return true;} return false;}


    private:
        int num_ice_like_atoms_;
        int num_ice_like_atoms_before_correction_;

        // whether we are performing surface correction or not 
        bool surface_correction_=false;

        // whether or not we are finding true ice using the InstantaneousInterface algorithm
        bool findtrueice_=false;
        bool findtrueclathrate_=false;
        Real n_, sigma_, isoval_;
        INT3 nL_;
        bool pbcMesh_=true;
        densityptr density_;
        Real3 volume_={{0,0,0}};
        bool cutMesh_=false;
        Real3 Ray_;

        // keep track of the water indices in the system 
        std::vector<int> water_indices_;
        std::vector<int> ice_like_indices_;
        std::vector<int> is_ice_like_;
        std::vector<int> is_clathrate_like_;
        std::vector<int> clathrate_like_indices_;
        Real surface_cutoff_=0.6;
        Real surface_cutoff_sq_;
        Real ice_cutoff_=0.55;
        Real ice_cutoff_sq_;
        int surface_threshold_=1;
        int ice_threshold_=8;

        // cell list related things 
        cellptr cell_;
        cellptr cell_surface_;
        cellptr cell_ice_corr_;
        cellptr cell_surface_corr_;
        std::vector<std::vector<int>> water_cell_list_;
        std::vector<std::vector<std::vector<int>>> surface_cell_list_;
        
        Real solvation_shell_r_=0.35, solvation_shell_r_squared_;
        Real surface_shell_r_=0.45, surface_shell_r_squared_;

        int harmonics_degree_=3, num_m_;
        std::vector<int> m_vec_;

        // atom group names 
        std::string atomgroup_name_;
        std::vector<std::string> surface_atomgroups_;

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

        std::vector<std::vector<int>> Ice_Indices_;

        // probe volume related things 
        std::vector<int> IsInsideProbeVolume_;

        // added ice index
        std::vector<int> addedIce_;

        // types 
        std::vector<int> types_;
        std::array<int, 7> type_nums_;
};
