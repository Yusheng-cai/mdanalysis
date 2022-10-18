#pragma once

#include "Calculation.h"
#include "Lattice.hpp"
#include "CellGrid.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"
#include "marching_cubes.hpp"
#include "parallel/OpenMP_buffer.h"
#include "Bin.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

class QtensorLattice : public Calculation{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using Matrix  = CommonTypes::Matrix;
        using INT2    = CommonTypes::index2;
        using Binptr  = std::unique_ptr<Bin>;

        QtensorLattice(const CalculationInput& input);

        virtual void calculate() override;
        void calculateCoraseGrain();
        void calculateNeighborSearch();
        virtual void update() override;
        virtual void finishCalculate() override;

        // initialize lattice offsets
        void CalculateLatticeOffsets();

        // calculate corase grain functionj
        Real CalculateCoraseGrainFunction(Real& rsq);

        // calculate azimuthal qtensor
        void CalculateZenithalQtensor();

        // printing function
        void printDirector(std::string name);
        void printOrder(std::string name);
        void printQtensor(std::string name);
        void printVelocity(std::string name);
        void printReducedDirector(std::string name);
        void printReducedOrder(std::string name);
        void printIsoSurface(std::string name);
        void printZenithalOrder(std::string name);
        void printAzimuthalOrder(std::string name);

        // per iter printing function
        void printOrderPerIter(std::ofstream& ofs);
        void printDensityPerIter(std::ofstream& ofs);

    private:
        bool coarse_grain_=false;

        // define the cellgrid
        cellptr cell_;

        // define the lattice 
        Lattice<Real3> lattice_;
        Lattice<Matrix> lattice_Qtensor_;
        Lattice<Real3> lattice_director_;
        Lattice<Real> lattice_order_;
        Lattice<Real> lattice_num_atoms_;
        Lattice<Real> lattice_biaxiality_;

        // per iter
        Lattice<Matrix> lattice_Qtensor_Iter_;
        Lattice<Real> lattice_num_atoms_Iter_;
        Lattice<Real> lattice_Order_Iter_;

        //          coarse graining parameters          //
        INT3 lattice_shape_;
        Real3 dL_;
        Real sigma_, sigma2_, prefactor_, inv_factor_;
        int n_;
        std::vector<INT3> lattice_offsets_;



        //          non coarse-graining parameters      //
        Real cutoff_sq_;
        Real min_dist_sq_;
        Real cutoff_;
        Real min_dist_;

        // resname
        std::string resname_;
        std::string refResname_;
        std::vector<int> refCOMIndices_;
        bool reference_=false;
        Real3 refCOM_;

        // zenithal order
        std::vector<std::vector<Matrix>> Zenithal_Qtensor_;
        std::vector<std::vector<Real>> Zenithal_num_;
        std::vector<std::vector<Real>> Zenithal_Order_;

        // azimuthal order
        std::vector<std::vector<Matrix>> Azimuthal_Qtensor_;
        std::vector<std::vector<Real>> Azimuthal_num_;
        std::vector<std::vector<Real>> Azimuthal_Order_;

        // R index 
        INT3 Rindex_ = {{0,1,2}};

        Binptr Rbin_;
        int numrbin_;
        Binptr Ztbin_;
        Binptr Atbin_;
        int numtbin_;
        bool usePredefinedDir_=false;
        Real3 predefinedDir_;
        Real3 rotate_vec_={{0,0,1}};

        // uij
        std::vector<Real3> uij_;

        // whether we are reducing the lattice on 2 dimensions
        INT2 reduced_dimensions_={{0,1}};
        bool reduced_=false;
        std::vector<std::vector<Matrix>> reduced_Q_;
        std::vector<std::vector<Real>> reduced_num_;
        std::vector<std::vector<Real3>> reduced_director_;
        std::vector<std::vector<Real>> reduced_order_;

        // head and tail index 
        int headIndex_=1, tailIndex_=2;

        // are we doing marching cubes
        MarchingCubes mc_;
        bool performMC_=false;
        std::string MC_name_;
        Real isoval_;
        Mesh m_;
        bool pbc_=true;

        // global Qtensor 
        Matrix GlobalQtensor_;
        Real3 GlobalDirector_;
        Real GlobalOrder_;
};