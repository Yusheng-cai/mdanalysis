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
#include <vector>
#include <array>
#include <complex>

class TrueIceFilter : public Calculation{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using INT2    = CommonTypes::index2;
        using densityptr = std::unique_ptr<DensityField>;

        TrueIceFilter(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        Real getNumIceLikeAtoms() const {return num_ice_like_atoms_;}
        Real getNumIceLikeChillPlus() const {return num_ice_like_atoms_before_correction_;}

        void CorrectIceLikeAtomsBasedOnSurface();
        void InterfaceFiltering();

        void printIceLikeAtoms(std::ofstream& ofs);
        void printAddedIce(std::ofstream& ofs);

    private:
        std::string filename_;

        int num_ice_like_atoms_;
        int num_ice_like_atoms_before_correction_;

        // keep track of the water indices in the system 
        std::vector<int> water_indices_;
        std::vector<std::vector<int>> ice_like_indices_total_;
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

        // whether or not we are finding true ice using the InstantaneousInterface algorithm
        bool InterfaceFiltering_=false;
        Real n_, sigma_, isoval_;
        INT3 nL_;
        bool pbcMesh_=true;
        densityptr density_;
        Real3 volume_={{0,0,0}};
        bool cutMesh_=false;
        Real3 Ray_;

        // cell list related things 
        cellptr cell_ice_corr_;
        cellptr cell_surface_corr_;
        std::vector<std::vector<int>> water_cell_list_;
        std::vector<std::vector<std::vector<int>>> surface_cell_list_;

        // atom group names 
        std::string atomgroup_name_;
        std::vector<std::string> surface_atomgroups_;

        // probe volume related things 
        std::vector<int> IsInsideProbeVolume_;

        // added ice index
        std::vector<int> addedIce_;
};
