#pragma once

#include "Lattice.hpp"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "parallel/OpenMP_buffer.h"
#include "marching_cubes.hpp"
#include "tools/Constants.h"
#include "tools/CommonOperations.h"
#include "Calculation.h"

#include <vector>
#include <map>
#include <array>
#include <string>
#include <memory>
#include <iostream>
#include <sstream>
#include <functional>

struct DensityFieldInput{
    using INT3 = CommonTypes::index3;
    using Real = CommonTypes::Real;

    SimulationState& simstate;
    INT3 nL;
    Real sigma;
    Real n;
    Real isoval;
    bool pbc;
};

class DensityField{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT3 =std::array<int,3>;
        using INT2 = CommonTypes::index2;
        using Range= CommonTypes::Real2;
        using Meshptr = std::unique_ptr<Mesh>;
        using outputfunc = std::function<void(std::string)>;

        DensityField(const DensityFieldInput& input);

        void update();

        // initialize the mesh object
        void initializeMesh();

        // helper function that calculates the offset index 
        void CalcOffsetIndex();

        // helper function to calculate the lattice point
        INT3 CalculateLatticePoint(const Real3& pos);

        // helper function to calculate the lattice position
        Real3 CalculateLatticePosition(const INT3& index);

        // calculate the instantaneous concentration field 
        void CalculateInstantaneousField(const std::vector<Real3>& pos, Mesh& mesh);

        // inline function that coarse grains the concentration
        inline Real GaussianCoarseGrainFunction(Real distsq);

    protected:
        // initialize the lattice 
        Lattice<Real> density_;

        // set up the openmp buffer for density field
        OpenMP::OpenMP_buffer<Lattice<Real>> density_buffer_;

        // The <dx,dy,dz> of the system
        Real3 dL_, L_;

        // The <Nx, Ny, Nz> of the system
        INT3 nL_;

        // The simulation state that keeps track of atom groups and bounding box
        SimulationState& simstate_;

        // sigma used for gaussian smoothing 
        Real sigma_, sigmasq_, prefactor_;
        // we cut off n sigmas away
        Real n_, cutoff_;

        // offset indices 
        std::vector<INT3> offsetIndex_;

        // The isosurface value, usually in units of atom/nm3
        Real isoSurfaceVal_;

        MarchingCubes MarchingCubes_;

        // Mesh stuff 
        bool pbcmesh_;
};