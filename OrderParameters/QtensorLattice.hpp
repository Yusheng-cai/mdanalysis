#pragma once

#include "Calculation.h"
#include "Lattice.hpp"
#include "CellGrid.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

class QtensorLattice : public Calculation{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        using Matrix  = CommonTypes::Matrix;
        using INT2    = CommonTypes::index2;

        QtensorLattice(const CalculationInput& input);

        virtual void calculate() override;
        virtual void update() override;
        virtual void finishCalculate() override;

        // printing function
        void printDirector(std::string name);
        void printOrder(std::string name);
        void printQtensor(std::string name);
        void printVelocity(std::string name);
        void printReducedDirector(std::string name);
        void printReducedOrder(std::string name);

    private:
        // define the cellgrid
        cellptr cell_;

        // define the lattice 
        Lattice<Real3> lattice_;
        Lattice<Matrix> lattice_Qtensor_;
        Lattice<Real3> lattice_director_;
        Lattice<Real> lattice_order_;
        Lattice<Real> lattice_num_atoms_;
        Lattice<Real> lattice_biaxiality_;

        INT3 lattice_shape_;
        Real3 dL_;

        // resname
        std::string resname_;

        // the cutoff
        Real cutoff_;
        Real cutoff_sq_;
        Real min_dist_=0.5;
        Real min_dist_sq_;

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
};