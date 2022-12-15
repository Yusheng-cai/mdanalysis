#pragma once
#include <cstdio>

#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"
#include "tools/Algorithm.h"
#include "OrderParameters/Lattice.hpp"
#include "OrderParameters/CellGrid.h"
#include "MorsePotential.hpp"
#include "LennardJones.hpp"
#include "parallel/OpenMP_buffer.h"
#include "OrderParameters/LinAlgTools.h"
#include "mda_tools.hpp"

#include <memory>

class NanoparticleGeneration
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT3 = CommonTypes::index3;
        using LJptr= std::unique_ptr<LennardJones>;
        using Mptr = std::unique_ptr<MorsePotential>;
        using cellptr = std::unique_ptr<CellGrid>;

        NanoparticleGeneration(const std::vector<Real3>& particle_pos, const std::string& ligand_gro, const std::string& ligand_top, Real De=36.664, Real alpha=14.7, \
                                                        Real re=0.265, Real sigma=0.425, Real epsilon=1.661);

        void ProcessParticlePositions();

        void GeneratePotentialGrid();

        INT3 FindIJKOnGrid(const Real3& pos);
        Real3 GridPositionFromIJK(const INT3& ijk);
        void FixIndex(INT3& ijk);

        void CalculateDistance(const Real3& v1, const Real3& v2, Real3& dist, Real& distsq);

        void Generate();

        // initialize ligand information
        void initializeLigandInformation();

        // 
        void addLigandInformation();

        void CreateOffset();

        // function that adds a sulfur onto the grid, updates the energies in terms of LJ
        void addSulfur(int idx);

        // find low energy site near minimum
        int find_low_energy_site_nearmin();

        // getters
        const std::vector<Real3>& getSulfurPositions() const {return sulfur_positions_;}

        // energy minimization
        void MinimizeEnergy();

        // construct final structure with ligand
        void constructNPWithLigand();

        void writeGroFile(std::string name);

    private:
        Lattice<Real> potential_grid_;
        OpenMP::OpenMP_buffer<Lattice<Real>> potential_grid_buffer_;

        std::vector<Real3> particle_pos_;

        Real3 box_size_, box_min_, box_max_;
        Real box_step_=0.1;
        INT3 box_N_;

        // Parameters for potential
        Real De_, alpha_, re_, sigma_, epsilon_;

        // Instantiate the Morse and Lennard Jones potential
        Mptr Morse_;
        LJptr LJ_;

        // number of gold particles
        int numParticles_;

        // cutoff 
        Real cutoff_ = 1.0;
        Real cutoff_sq_ = 1.0;

        // offset 
        std::vector<INT3> offsets_;

        // sulfur positions
        std::vector<Real3> sulfur_positions_;

        // minimum energy as well as index
        int min_index_;
        Real min_energy_;

        // threshold energy
        Real e_thresh_ = -50; // kJ/mol

        // grofile c string
        const char* gro_c_string="%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n";

        // ligand gro and ligand top files
        std::string ligand_gro_, ligand_top_;

        // ligand_gro positions
        std::vector<Real3> ligand_positions_;
        std::vector<std::string> ligand_names_, ligand_resnames_;
};