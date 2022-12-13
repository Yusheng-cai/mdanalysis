#pragma once

#include "tools/CommonTypes.h"
#include "gemmi/cif.hpp"
#include "gemmi/numb.hpp"
#include "tools/Assert.h"
#include "tools/CommonOperations.h"
#include "tools/Constants.h"
#include "OrderParameters/LinAlgTools.h"
#include "tools/Algorithm.h"
#include "mda_tools.hpp"

#include <string>
#include <cmath>
#include <vector>

class Crystal
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;

        Crystal(std::string fname);

        void construct_matrix();

        void calculate();

        void writeOutput(const std::string& fname);


    private:
        std::string cif_fname_;

        // positions
        std::vector<Real3> cell_positions_;

        // a,b,c,alpha,gamma,omega
        Real a_,b_,c_,alpha_,beta_,gamma_;

        // construct the Matrix
        Matrix m_;

        // nx, ny, nz
        int nx_=2, ny_=2, nz_=1;

        // total cell positions
        std::vector<Real3> total_cell_pos_;
        std::vector<std::string> total_atom_types_;

        std::vector<std::string> unit_cell_atom_types_;
        std::vector<std::string> unique_atom_types_;
        std::vector<std::string> total_resname_;
};