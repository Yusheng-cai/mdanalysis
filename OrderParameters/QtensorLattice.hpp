#pragma once

#include "Calculation.h"
#include "Lattice.hpp"
#include "CellGrid.h"
#include "SimulationState.h"
#include "LinAlgTools.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

class QtensorLattice : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;

        QtensorLattice(const CalculationInput& input);

        virtual void calculate();
        virtual void update();

    private:
        // define the cellgrid
        cellptr cell_;

        // define the lattice 
        Lattice<Real3> lattice_;
        INT3 lattice_shape_;
        int nx, ny, nz;
        Real3 dL_;

        // resname
        std::string resname_;

        // the cutoff
        Real cutoff_;

        // uij
        std::vector<Real3> uij_;

        // head and tail index 
        int headIndex_=1, tailIndex_=2;
};